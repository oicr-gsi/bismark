package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class BismarkMergeSortDedup extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outputDir;
    private String sample;
    private String sampleName;

    // IO Data
    private String bamFiles;
    private Map<String, Set<String>> inputFilesByGroup;
    private List<String> outputIdentifiers;
    private final Map<String, Set<String>> provisionedFilesByGroup = new HashMap<>();

    // Output check
    private boolean isFolder = true;

    //Scripts 
    private String bismarkDedup;

    //Tools
    private String samtools;

    //Memory allocation
    private Integer bismarkMem;

    // bismark params
    private Integer threads;
    private Integer seedLength;
    private Integer cores;
    private Integer maxMismatch;

    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(BismarkMergeSortDedup.class.getName());
    private String queue;
    private Map<String, SqwFile> tempFiles;

    // provision
    private final static String BAM_METATYPE = "application/bam";
    private final static String BAI_METATYPE = "application/bam-index";
    private final static String TXT_METATYPE = "text/plain";

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // output
            //output identifiers
            outputIdentifiers = Arrays.asList(StringUtils.split(getProperty("output_identifiers"), ";"));
            if (hasDuplicates(outputIdentifiers)) {
                error("Duplicates detected in output_identifiers");
            }

            // input samples 
            bamFiles = getProperty("input_bam_files");
            inputFilesByGroup = getMapOfSets(outputIdentifiers, Arrays.asList(StringUtils.split(getProperty("input_bam_files"), ";")), ",");

            // output filename prefix ; if "output_filename_prefix exists parse from here
            sampleName = getProperty("output_filename_prefix");
            //samtools
            samtools = getProperty("samtools");

            //bismark utils
            bismarkDedup = getProperty("bismark_deduplicate");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            //memory
            bismarkMem = Integer.parseInt(getProperty("bismark_mem"));

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        int id = 0;
        for (Map.Entry<String, Set<String>> e : inputFilesByGroup.entrySet()) {
            String outputGroup = e.getKey();
            Set<String> provisionedPaths = new HashSet<>();
            for (String filePath : e.getValue()) {
                if (!"bam".equals(FilenameUtils.getExtension(filePath))) {
                    error("Unsupported input file: " + filePath);
                }
                SqwFile bam = this.createFile("file_in_" + id++);
                bam.setSourcePath(filePath);
                bam.setType(BAM_METATYPE);
                bam.setIsInput(true);
                provisionedPaths.add(bam.getProvisionedPath());
            }
            provisionedFilesByGroup.put(outputGroup, provisionedPaths);
        }

        return this.getFiles(); // I have an array of input files
    }

    @Override
    public void buildWorkflow() {
        /**
         * This workflow reads in the fastq files generated through Methyl-Seq
         * experiments Runs bismark to convert C-> T and G-> A (in reverse
         * strand) in the fastq files (to identify methylation sites) Aligns the
         * intermediate fastqfiles to reference sequence (using bowtie2 and
         * converts output to sam format using samtools) Runs bismark
         * methylation extractor to identify methylation sites
         */

        this.outputDir = this.dataDir + "output/";

        //setup operations performed on output group map
        Map<String, String> outputNameByOutputGroup = new HashMap<>();
        for (String outputIdentifier : outputIdentifiers) {
            outputNameByOutputGroup.put(outputIdentifier, outputIdentifier + ".");
        }

        //used to collect input files (and associated job) to the next step of the workflow - partitioned by the output group
        Multimap<String, Pair<String, Job>> inputFileAndJobToNextStepByOutputGroup = null;

        //for each output group: merge, sort, filter, dedup, index
        Multimap<String, Pair<String, Job>> mergedBamsByGroup = HashMultimap.create();
        for (String outputGroup : outputIdentifiers) {

            Job parentJob;

            if (inputFilesByGroup.get(outputGroup).size() > 1) {
                this.sampleName = outputNameByOutputGroup.get(outputGroup);
//                String mergeSortedBam = this.sampleName + "sorted.";
                outputNameByOutputGroup.put(outputGroup, this.sampleName);

                String[] filemap = provisionedFilesByGroup.get(outputGroup).toArray(new String[0]);
                ArrayList<String> inputFiles = new ArrayList<String>();
                for (String b : filemap) {
                    inputFiles.add(b);
                }
                Job samtoolsMergeSortBam = samtoolsMerge(inputFiles);
                parentJob = samtoolsMergeSortBam;

                // Expected intermediate outputs 
                String mergedBam = this.outputDir + "/" + this.sampleName + ".bam";
                String mergeSortedBam = this.outputDir + "/" + this.sampleName + "sorted.bam";

                // sort merged bam file
                Job jobSamtoolsSort = samtoolsSort(mergedBam, mergeSortedBam);
                jobSamtoolsSort.addParent(parentJob);
                parentJob = jobSamtoolsSort;

                // Expected dedup bam filename 
                String dedupBam = this.outputDir + "/" + this.sampleName + ".sorted.deduplicated.bam";

//                deduplicate bam
                Job deduplicateBam = deduplicateBismarkBam(mergeSortedBam);
                deduplicateBam.addParent(parentJob);
                parentJob = deduplicateBam;

                // index dedup file
                Job indexBam = indexBam(dedupBam);
                indexBam.addParent(parentJob);
                parentJob = indexBam;

                // provision bam outputs
                SqwFile bamFile = createOutputFile(dedupBam, BAM_METATYPE, this.manualOutput);
                bamFile.getAnnotations().put("deduplicated bam file ", "bismark_dedup");
                deduplicateBam.addFile(bamFile);

                // provision bai outputs
                SqwFile baiFile = createOutputFile(dedupBam + ".bai", BAI_METATYPE, this.manualOutput);
                bamFile.getAnnotations().put("index file of deduplicated bam ", "bismark_dedup");
                indexBam.addFile(baiFile);
                
                // provision dedup report file
                SqwFile reportFile = createOutputFile(dedupBam.replace(".bam", "_report.txt"), TXT_METATYPE, this.manualOutput);
                reportFile.getAnnotations().put("report file", "bismark_dedup");
                deduplicateBam.addFile(reportFile);
            }
        }

    }

    private Job samtoolsMerge(ArrayList<String> bamFiles) {
        // given an array of bam files, merge them
        String inBams = bamFiles.toString();
        Job jobMergeBam = getWorkflow().createBashJob("merge_bams");
        Command cmd = jobMergeBam.getCommand();
        cmd.addArgument("export PATH=" + this.samtools + ":$PATH;");
        cmd.addArgument(this.samtools + "/samtools merge ");
        cmd.addArgument(inBams); // output merge bam
        jobMergeBam.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobMergeBam.setQueue(getOptionalProperty("queue", ""));
        return jobMergeBam;
    }

    private Job deduplicateBismarkBam(String sortedBam) {
        Job jobDedupBam = getWorkflow().createBashJob("bismark_deduplicate");
        Command cmd = jobDedupBam.getCommand();
        cmd.addArgument(bismarkDedup);
        cmd.addArgument("-p --bam " + sortedBam);
        jobDedupBam.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobDedupBam.setQueue(getOptionalProperty("queue", ""));
        return jobDedupBam;
    }

    private Job samtoolsSort(String mergeBam, String mergeSortedBam) {
        // sort bamFile
        Job jobSortBam = getWorkflow().createBashJob("sort_bam");
        Command cmd = jobSortBam.getCommand();
        cmd.addArgument("export PATH=" + this.samtools + ":$PATH;");
        cmd.addArgument(this.samtools + "/samtools sort");
        cmd.addArgument("-o " + mergeSortedBam);
        cmd.addArgument("-O bam");
        cmd.addArgument("-T " + this.tmpDir + "/" + this.sampleName + "_tmp");
        cmd.addArgument(mergeBam);
        jobSortBam.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobSortBam.setQueue(getOptionalProperty("queue", ""));
        return jobSortBam;
    }

    private Job indexBam(String inputBam) {
        // bismark methylation extractor
        Job jobToBai = getWorkflow().createBashJob("index_bam");
        Command cmd = jobToBai.getCommand();
        cmd.addArgument("export PATH=" + this.samtools + ":$PATH;");
        cmd.addArgument(this.samtools + "/samtools index " + inputBam);
        jobToBai.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobToBai.setQueue(getOptionalProperty("queue", ""));
        return jobToBai;
    }

    private Map<String, Set<String>> getMapOfSets(List<String> keys, List<String> values, String valuesDelimiter) {
        if (keys.isEmpty() || values.isEmpty() || keys.size() != values.size()) {
            error("key size (" + keys.size() + ") != values size (" + values.size() + ")");
        }
        Map<String, Set<String>> map = new HashMap<>();
        for (int i = 0; i < keys.size(); i++) {
            List<String> valsList = Arrays.asList(values.get(i).split(valuesDelimiter));
            Set<String> vals = new LinkedHashSet<>(valsList);
            if (map.put(keys.get(i), vals) != null) {
                error("Duplicate key detected = [" + keys.get(i));
            }
        }
        return map;
    }

    private <T> boolean hasDuplicates(Collection<T> c) {
        Set<T> s = new HashSet<>(c);
        return c.size() != s.size();
    }

    private void error(String msg) {
        Log.error(msg);
        setWorkflowInvalid();
    }

}
