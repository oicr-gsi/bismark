package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.io.File;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

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
public class WorkflowClient extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;

    // Input Data
    private String fastqFiles;
    private String genomeFolder;

    // Output check
    private boolean isFolder = true;

    //Scripts 
    private String bismark;
    private String bismarkMethylationExtractor;

    //Tools
    private String samtools;
    private String bowtie;

    //Memory allocation
    private Integer bismarkMem;

    // bismark params
    private Integer threads;
    private Integer seedLength;
    private Integer cores;
    private Integer maxMismatch;


    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(WorkflowClient.class.getName());
    private String queue;
    private Map<String, SqwFile> tempFiles;

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            fastqFiles = getProperty("fastq_files");

            //samtools
            samtools = getProperty("samtools");

            //bowtie
            bowtie = getProperty("bowtie");

            // genome folder
            genomeFolder = getProperty("genome_folder");

            //bismark utils
            bismark = getProperty("bismark");
            bismarkMethylationExtractor = getProperty("bismark_methylation_extractor");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            //memory
            bismarkMem = Integer.parseInt(getProperty("bismark_mem"));

            // params
            maxMismatch = Integer.parseInt(getProperty("max_mismatch_allowed"));
            seedLength = Integer.parseInt(getProperty("seed_length"));
            threads = Integer.parseInt(getProperty("no_of_threads"));
            cores = Integer.parseInt(getProperty("no_of_multiprocessing_cores"));


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
    public void buildWorkflow() {
        /**
         * This workflow reads in the fastq files generated through Methyl-Seq
         * experiments Runs bismark to convert C-> T and G-> A (in reverse
         * strand) in the fastq files (to identify methylation sites) Aligns the
         * intermediate fastqfiles to reference sequence (using bowtie2 and
         * converts output to sam format using samtools) Runs bismark
         * methylation extractor to identify methylation sites
         */
        Job parentJob = null;
        String[] fqFile = fastqFiles.split(",");
        String outputDir = "output/";
        String[] pathsplit = fqFile[0].split("/");
        Integer n = pathsplit.length;
        String name = pathsplit[n - 1];
        String[] names = name.split("\\.");
        String sample = names[0];
        String expectedOutputSam = sample + "_pe.sam";

        Job bismark = runBismarkPreprocess(fqFile, sample, outputDir);
        parentJob = bismark;

        // check if the out directory contains the expectedOutputSam file -- can go in decider
        Job jobBismarkMethylationExtractor = jobMethylationExtractor(expectedOutputSam, outputDir);
        jobBismarkMethylationExtractor.addParent(parentJob);
        parentJob = jobBismarkMethylationExtractor;
    }

    // create Job function for the bismark alignment
    private Job runBismarkPreprocess(String[] fastqFiles, String sample, String outDir) {
        String fastq1Path = fastqFiles[0];
        String fastq2Path = fastqFiles[1];
        Job jobBismark = getWorkflow().createBashJob("bismark");
        Command command = jobBismark.getCommand();
        command.addArgument(bismark);
        command.addArgument("--path_to_bowtie " + bowtie);
        command.addArgument("--sam");
        command.addArgument("-n " + maxMismatch.toString());
        command.addArgument("-l " + seedLength.toString());
        command.addArgument("-p " + threads.toString()); //l -> seed length; -n --> max no. of mismatched permitted in the seed; p-number of threads to parallelize the job
        command.addArgument(this.genomeFolder);
        command.addArgument("-b " + sample);
        command.addArgument("-o " + outDir);
        command.addArgument("--temp_dir " + tmpDir);
        command.addArgument("--gzip");
        command.addArgument("-1 " + fastq1Path);
        command.addArgument("-2 " + fastq2Path);
        jobBismark.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobBismark.setQueue(getOptionalProperty("queue", ""));
        return jobBismark;
    }

    private Job jobMethylationExtractor(String samFile, String outDir) {
        // bismark methylation extractor
        Job jobMethExtractor = getWorkflow().createBashJob("methylation_extractor");
        Command cmd = jobMethExtractor.getCommand();
        cmd.addArgument(bismarkMethylationExtractor);
        cmd.addArgument("-p"); //paired end
        cmd.addArgument("--comprehensive");
        cmd.addArgument("--merge_non_CpG");
        cmd.addArgument("-o " + outDir);
        cmd.addArgument("--genome_folder " + genomeFolder);
        cmd.addArgument("--samtools_path " + samtools);
        cmd.addArgument("--gzip");
        cmd.addArgument("--multicore " + cores.toString()); //no. of parallel instances bismark to run
        cmd.addArgument("--bedGraph");
        cmd.addArgument("--cytosine_report " + outDir + '/' + samFile);
        jobMethExtractor.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobMethExtractor.setQueue(getOptionalProperty("queue", ""));
        return jobMethExtractor;
    }

}
