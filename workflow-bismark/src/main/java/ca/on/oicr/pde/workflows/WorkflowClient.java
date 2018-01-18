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
    private String expectedOutputBam;
    private String outputDir;
    private String sample;
    private String sampleName;

    // Input Data
    private String r1FastqFile;
    private String r2FastqFile;
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

    // provision
    private final static String FASTQ_METATYPE = "chemical/seq-na-fastq-gzip";
    private final static String BAM_METATYPE = "application/bam";
    private final static String BAI_METATYPE = "application/bam-index";
    private final static String TXT_GZ_METATYPE = "application/txt-gz";
    private final static String TXT_METATYPE = "text/plain";

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            r1FastqFile = getProperty("r1_fastq_file");
            r2FastqFile = getProperty("r2_fastq_file");
            
            // output filename prefix ; if "output_filename_prefix exists parse from here
            sampleName = getProperty("output_filename_prefix");
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
    public Map<String, SqwFile> setupFiles() {
        String r1FqFile = this.r1FastqFile;
        SqwFile file0 = this.createFile("R1");
        file0.setSourcePath(r1FqFile);
        file0.setType(FASTQ_METATYPE);
        file0.setIsInput(true);
        String r2FqFile = this.r2FastqFile;
        SqwFile file1 = this.createFile("R2");
        file1.setSourcePath(r2FqFile);
        file1.setType(FASTQ_METATYPE);
        file1.setIsInput(true);
        return this.getFiles();
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
        String r1FqFile = this.r1FastqFile;
        this.outputDir = this.dataDir + "output/";
        if (this.sampleName != null){
            String[] pathsplit = r1FqFile.split("/");
            Integer n = pathsplit.length;
            String name = pathsplit[n - 1];
            String[] names = name.split("\\.");
            this.sample = names[0];
        } else {
            this.sample = this.sampleName;
        }
//        this.sample = this.sampleName;
        this.expectedOutputBam = this.sample + "_pe.sam";

        Job bismark = runBismarkPreprocess();
        parentJob = bismark;

        // check if the out directory contains the expectedOutputBam file -- can go in decider
        Job jobBismarkMethylationExtractor = jobMethylationExtractor();
        jobBismarkMethylationExtractor.addParent(parentJob);
        parentJob = jobBismarkMethylationExtractor;

        // provision output files from BisMarkMethylationExtractor Job
        SqwFile cpgFile = createOutputFile(this.outputDir + "CpG_context_" + this.sample + "_pe.txt.gz", TXT_GZ_METATYPE, this.manualOutput);
        cpgFile.getAnnotations().put("CpG context file", "Bismark_methylation_extractor");
        parentJob.addFile(cpgFile);

        SqwFile nonCpGFile = createOutputFile(this.outputDir + "Non_CpG_context_" + this.sample + "_pe.txt.gz", TXT_GZ_METATYPE, this.manualOutput);
        nonCpGFile.getAnnotations().put("Non-CpG context file", "Bismark_methylation_extractor");
        parentJob.addFile(nonCpGFile);

        SqwFile splittingReportFile = createOutputFile(this.outputDir + this.sample + "_pe_splitting_report.txt", TXT_METATYPE, this.manualOutput);
        splittingReportFile.getAnnotations().put("splitting report", "Bismark_methylation_extractor");
        parentJob.addFile(splittingReportFile);

        SqwFile mBiasFile = createOutputFile(this.outputDir + this.sample + "_pe.M-bias.txt", TXT_METATYPE, this.manualOutput);
        mBiasFile.getAnnotations().put("M-bias", "Bismark_methylation_extractor");
        parentJob.addFile(mBiasFile);

        SqwFile bedGraphFile = createOutputFile(this.outputDir + this.sample + "_pe.bedGraph.gz", TXT_GZ_METATYPE, this.manualOutput);
        bedGraphFile.getAnnotations().put("bedGraph", "Bismark_methylation_extractor");
        parentJob.addFile(bedGraphFile);

        SqwFile bismarkCovFile = createOutputFile(this.outputDir + this.sample + "_pe.bismark.cov.gz", TXT_GZ_METATYPE, this.manualOutput);
        bismarkCovFile.getAnnotations().put("coverage report", "Bismark_methylation_extractor");
        parentJob.addFile(bismarkCovFile);

        SqwFile cpgReport = createOutputFile(outputDir + this.sample + "_pe.CpG_report.txt.gz", TXT_GZ_METATYPE, this.manualOutput);
        cpgReport.getAnnotations().put("CpG report", "Bismark_methylation_extractor");
        parentJob.addFile(cpgReport);

//        // convert Sam to Bam
//        Job jobCovertSam = jobSamToBam();
//        jobCovertSam.addParent(bismark);
//        //parentJob = jobCovertSam;
//
//        Job jobSortBamFile = jobSortBam();
//        jobSortBamFile.addParent(jobCovertSam);
        //parentJob = jobSortBamFile;
        //String sortedBamFile = this.expectedOutputBam.replace(".sam", ".sorted.bam");
        Job jobIndexSortedBamFile = jobIndexBam();
        jobIndexSortedBamFile.addParent(parentJob);
        //parentJob = jobIndexBamFile;

        // provision bam outputs
        SqwFile bamFile = createOutputFile(this.outputDir + this.expectedOutputBam.replace(".sam", ".sorted.bam"), BAM_METATYPE, this.manualOutput);
        bamFile.getAnnotations().put("segment data from the tool ", "Sequenza ");
        parentJob.addFile(bamFile);

        SqwFile baiFile = createOutputFile(this.outputDir + this.expectedOutputBam.replace(".sam", ".sorted.bam.bai"), BAI_METATYPE, this.manualOutput);
        bamFile.getAnnotations().put("segment data from the tool ", "Sequenza ");
        jobIndexSortedBamFile.addFile(baiFile);

    }

    // create Job function for the bismark alignment
    private Job runBismarkPreprocess() {
        //String fastq1Path = fastqFiles[0];
        //String fastq2Path = fastqFiles[1];
        Job jobBismark = getWorkflow().createBashJob("bismark");
        Command command = jobBismark.getCommand();
        command.addArgument(bismark);
        command.addArgument("--path_to_bowtie " + bowtie);
        command.addArgument("--bam");
        command.addArgument("--rg_tag"); // to add read group to bam
        command.addArgument("-n " + maxMismatch.toString());
        command.addArgument("-l " + seedLength.toString());
        command.addArgument("-p " + threads.toString()); //l -> seed length; -n --> max no. of mismatched permitted in the seed; p-number of threads to parallelize the job
        command.addArgument(this.genomeFolder);
        command.addArgument("-b " + this.sample);
        command.addArgument("-o " + this.outputDir);
        command.addArgument("--temp_dir " + tmpDir);
        command.addArgument("--gzip");
        command.addArgument("-1 " + getFiles().get("R1").getProvisionedPath());
        command.addArgument("-2 " + getFiles().get("R2").getProvisionedPath());
        jobBismark.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobBismark.setQueue(getOptionalProperty("queue", ""));
        return jobBismark;
    }

    private Job jobMethylationExtractor() {
        // bismark methylation extractor
        Job jobMethExtractor = getWorkflow().createBashJob("methylation_extractor");
        Command cmd = jobMethExtractor.getCommand();
        cmd.addArgument(bismarkMethylationExtractor);
        cmd.addArgument("-p"); //paired end
        cmd.addArgument("--comprehensive");
        cmd.addArgument("--merge_non_CpG");
        cmd.addArgument("-o " + this.outputDir);
        cmd.addArgument("--genome_folder " + genomeFolder);
        cmd.addArgument("--samtools_path " + samtools);
        cmd.addArgument("--gzip");
        cmd.addArgument("--multicore " + cores.toString()); //no. of parallel instances bismark to run
        cmd.addArgument("--bedGraph");
        cmd.addArgument("--cytosine_report " + this.outputDir + '/' + this.expectedOutputBam);
        jobMethExtractor.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobMethExtractor.setQueue(getOptionalProperty("queue", ""));
        return jobMethExtractor;
    }

//    private Job jobSamToBam() {
//        // bismark methylation extractor
//        String bamFile = this.expectedOutputBam.replace(".sam", ".bam");
//        Job jobToBam = getWorkflow().createBashJob("sam_to_bam");
//        Command cmd = jobToBam.getCommand();
//        cmd.addArgument(this.samtools + " view -bS " + this.outputDir + this.expectedOutputBam + " > " + this.outputDir + bamFile);
//        jobToBam.setMaxMemory(Integer.toString(bismarkMem * 1024));
//        jobToBam.setQueue(getOptionalProperty("queue", ""));
//        return jobToBam;
//    }
//
//    private Job jobSortBam() {
//        // bismark methylation extractor
//        String bamFile = this.expectedOutputBam.replace(".sam", ".bam");
//        Job jobSortBam = getWorkflow().createBashJob("sort_bam");
//        Command cmd = jobSortBam.getCommand();
//        cmd.addArgument(this.samtools + " sort " + this.outputDir + bamFile);
//        cmd.addArgument(this.outputDir + bamFile.replace(".bam", ".sorted"));
//        jobSortBam.setMaxMemory(Integer.toString(bismarkMem * 1024));
//        jobSortBam.setQueue(getOptionalProperty("queue", ""));
//        return jobSortBam;
//    }

    private Job jobIndexBam() {
        // bismark methylation extractor
        Job jobToBai = getWorkflow().createBashJob("index_bam");
        Command cmd = jobToBai.getCommand();
        cmd.addArgument(this.samtools + " index " + this.outputDir + 
                            this.expectedOutputBam); //.replace(".sam", ".sorted.bam"));
        jobToBai.setMaxMemory(Integer.toString(bismarkMem * 1024));
        jobToBai.setQueue(getOptionalProperty("queue", ""));
        return jobToBai;
    }

}
