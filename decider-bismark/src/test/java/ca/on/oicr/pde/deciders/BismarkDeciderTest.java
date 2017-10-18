//package ca.on.oicr.pde.deciders;
//
//import ca.on.oicr.gsi.provenance.DefaultProvenanceClient;
//import ca.on.oicr.gsi.provenance.ExtendedProvenanceClient;
//import ca.on.oicr.gsi.provenance.FileProvenanceFilter;
//import ca.on.oicr.gsi.provenance.ProvenanceClient;
//import ca.on.oicr.gsi.provenance.SeqwareMetadataAnalysisProvenanceProvider;
//import ca.on.oicr.gsi.provenance.SeqwareMetadataLimsMetadataProvenanceProvider;
//import ca.on.oicr.gsi.provenance.model.SampleProvenance;
//import ca.on.oicr.pde.client.SeqwareClient;
//import ca.on.oicr.pde.dao.executor.SeqwareExecutor;
//import ca.on.oicr.pde.testing.metadata.RegressionTestStudy;
//import ca.on.oicr.pde.model.SeqwareObject;
//import ca.on.oicr.pde.reports.WorkflowReport;
//import ca.on.oicr.pde.testing.metadata.RegressionTestStudy.SeqwareObjects;
//import ca.on.oicr.pde.utilities.Helpers;
//import ca.on.oicr.pde.testing.metadata.SeqwareTestEnvironment;
//import com.google.common.collect.ImmutableSet;
//import com.google.common.collect.Iterables;
//import com.google.common.collect.Sets;
//import java.io.File;
//import java.io.IOException;
//import java.util.Arrays;
//import java.util.Collection;
//import java.util.Collections;
//import java.util.HashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//import net.sourceforge.seqware.common.metadata.Metadata;
//import net.sourceforge.seqware.common.model.IUS;
//import net.sourceforge.seqware.common.model.Workflow;
//import net.sourceforge.seqware.common.module.FileMetadata;
//import org.apache.commons.lang3.ArrayUtils;
//import org.testng.Assert;
//import static org.testng.Assert.assertNotNull;
//import static org.testng.Assert.assertTrue;
//import org.testng.annotations.AfterClass;
//import org.testng.annotations.BeforeClass;
//import org.testng.annotations.BeforeSuite;
//import org.testng.annotations.Test;
//
///**
// *
// * @author mlaszloffy
// */
//public class BismarkDeciderTest {
//
//    private RegressionTestStudy r;
//    private File bundledWorkflow;
//    private ProvenanceClient provenanceClient;
//    private ExtendedProvenanceClient extendedProvenanceClient;
//    private SeqwareClient seqwareClient;
//    private SeqwareExecutor seqwareExecutor;
//    private SeqwareObjects seqwareObjects;
//    private SeqwareTestEnvironment testEnv;
//    private Metadata metadata;
//    private Map<String, String> config;
//
//    @BeforeSuite
//    public void setupMetadb() {
//
//        //get the workflow bundle associated with the decider
//        bundledWorkflow = Helpers.getBundledWorkflow();
//        Assert.assertNotNull(bundledWorkflow, "Unable to locate the workflow bundle.");
//        Assert.assertTrue(bundledWorkflow.exists(), "The workflow bundle [" + bundledWorkflow + "] does not exist.");
//
//        //get the database settings
//        String dbHost = System.getProperty("dbHost");
//        String dbPort = System.getProperty("dbPort");
//        String dbUser = System.getProperty("dbUser");
//        String dbPassword = System.getProperty("dbPassword");
//        assertNotNull(dbHost, "Set dbHost to a testing Postgres database host name");
//        assertNotNull(dbPort, "Set dbPort to a testing Postgres database port");
//        assertNotNull(dbUser, "Set dbUser to a testing Postgres database user name");
//
//        //get the seqware webservice war
//        String seqwareWarPath = System.getProperty("seqwareWar");
//        assertNotNull(seqwareWarPath, "seqwareWar is not set.");
//        File seqwareWar = new File(seqwareWarPath);
//        assertTrue(seqwareWar.exists(), "seqware was is not accessible.");
//
//        testEnv = new SeqwareTestEnvironment(dbHost, dbPort, dbUser, dbPassword, seqwareWar);
//
//        //get the regression test study and PDE's service objects
//        r = new RegressionTestStudy(testEnv.getSeqwareLimsClient());
//        seqwareClient = testEnv.getSeqwareClient();
//        seqwareExecutor = testEnv.getSeqwareExecutor();
//        seqwareObjects = r.getSeqwareObjects();
//
//        metadata = testEnv.getMetadata();
//        config = testEnv.getSeqwareConfig();
//
//        SeqwareMetadataLimsMetadataProvenanceProvider seqwareMetadataLimsMetadataProvenanceProvider = new SeqwareMetadataLimsMetadataProvenanceProvider(testEnv.getMetadata());
//        SeqwareMetadataAnalysisProvenanceProvider seqwareMetadataAnalysisProvenanceProvider = new SeqwareMetadataAnalysisProvenanceProvider(testEnv.getMetadata());
//        DefaultProvenanceClient dpc = new DefaultProvenanceClient();
//        dpc.registerAnalysisProvenanceProvider("seqware", seqwareMetadataAnalysisProvenanceProvider);
//        dpc.registerSampleProvenanceProvider("seqware", seqwareMetadataLimsMetadataProvenanceProvider);
//        dpc.registerLaneProvenanceProvider("seqware", seqwareMetadataLimsMetadataProvenanceProvider);
//        provenanceClient = dpc;
//        extendedProvenanceClient = dpc;
//    }
//
//    @BeforeClass
//    public void addFiles() {
//        //Install base caller to generate fastqs
//        Workflow upstreamWorkflow = seqwareClient.createWorkflow("CASAVA", "0.0", "");
//
//        //get all ius objects and link two fastqs
//        for (int i = 1; i <= 22; i++) {
//            SeqwareObject obj = seqwareObjects.get("IUS" + i);
//
//            Map<FileProvenanceFilter, Set<String>> filters = new HashMap<>();
//            filters.put(FileProvenanceFilter.sample, ImmutableSet.of(obj.getSwAccession().toString()));
//            Collection<SampleProvenance> sps = provenanceClient.getSampleProvenance(filters);
//            SampleProvenance sp = Iterables.getOnlyElement(sps);
//            IUS ius = seqwareClient.addLims("seqware", sp.getSampleProvenanceId(), sp.getVersion(), sp.getLastModified());
//
//            FileMetadata file1 = new FileMetadata();
//            file1.setDescription("description");
//            file1.setMd5sum("md5sum");
//            file1.setFilePath("/tmp/" + obj.getSwAccession() + ".R1.fastq.gz");
//            file1.setMetaType("chemical/seq-na-fastq-gzip");
//            file1.setType("type?");
//            file1.setSize(1L);
//
//            FileMetadata file2 = new FileMetadata();
//            file2.setDescription("description");
//            file2.setMd5sum("md5sum");
//            file2.setFilePath("/tmp/" + obj.getSwAccession() + ".R2.fastq.gz");
//            file2.setMetaType("chemical/seq-na-fastq-gzip");
//            file2.setType("type?");
//            file2.setSize(1L);
//
//            seqwareClient.createWorkflowRun(upstreamWorkflow, ImmutableSet.of(ius), Collections.EMPTY_LIST, Arrays.asList(file1, file2));
//        }
//    }
//
//    @Test
//    public void basicScheduleTest() throws IOException {
//        Workflow bismarkWorkflow = seqwareExecutor.installWorkflow(Helpers.getBundledWorkflow());
//
//        run(bismarkWorkflow, "--all");
//        WorkflowReport report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report.getMaxInputFiles().intValue(), 2);
//        Assert.assertEquals(report.getMinInputFiles().intValue(), 2);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 16);
//        Assert.assertEquals(report.getFileMetaTypes(), Sets.newHashSet("chemical/seq-na-fastq-gzip"));
//
//        //all fastqs have been scheduled for processing, no new workflow runs should be schedulable
//        run(bismarkWorkflow, "--all");
//        WorkflowReport report2 = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report2.getWorkflowRunCount().intValue(), 16);
//        Assert.assertEquals(report, report2);
//    }
//
//    @Test
//    public void launchMaxTest() throws IOException {
//        Workflow bismarkWorkflow = seqwareExecutor.installWorkflow(Helpers.getBundledWorkflow());
//
//        WorkflowReport report;
//
//        //schedule 4 runs
//        run(bismarkWorkflow, "--all", "--launch-max", "4");
//        report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 4);
//
//        //schedule 4 more runs
//        run(bismarkWorkflow, "--all", "--launch-max", "4");
//        report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 8);
//
//        //shouldn't schedule anything
//        run(bismarkWorkflow, "--all", "--launch-max", "0");
//        report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 8);
//
//        //schedule the rest
//        run(bismarkWorkflow, "--all", "--launch-max", "20");
//        report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 16);
//    }
//
//    @Test
//    public void missingRootSampleFilterTest() throws IOException {
//        Workflow bismarkWorkflow = seqwareExecutor.installWorkflow(Helpers.getBundledWorkflow());
//
//        run(bismarkWorkflow, "--root-sample-name", "does_not_exist");
//
//        WorkflowReport report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 0);
//    }
//
//    @Test
//    public void missingSequencerRunFilterTest() throws IOException {
//        Workflow fastqcWorkflow = seqwareExecutor.installWorkflow(Helpers.getBundledWorkflow());
//
//        run(fastqcWorkflow, "--sequencer-run-name", "does_not_exist");
//
//        WorkflowReport report = WorkflowReport.generateReport(seqwareClient, provenanceClient, fastqcWorkflow);
//        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 0);
//    }
//
////    @Test
////    public void missingStudyFilterTest() throws IOException {
////        Workflow bismarkWorkflow = seqwareExecutor.installWorkflow(Helpers.getBundledWorkflow());
////
////        run(bismarkWorkflow, "--study-name", "does_not_exist");
////
////        WorkflowReport report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
////        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 0);
////    }
//
////    @Test
////    public void missingSampleFilterTest() throws IOException {
////        Workflow bismarkWorkflow = seqwareExecutor.installWorkflow(Helpers.getBundledWorkflow());
////
////        run(bismarkWorkflow, "--sample-name", "does_not_exist");
////
////        WorkflowReport report = WorkflowReport.generateReport(seqwareClient, provenanceClient, bismarkWorkflow);
////        Assert.assertEquals(report.getWorkflowRunCount().intValue(), 0);
////    }
//
////    private void run(Workflow workflow, String... extraDeciderParams) throws IOException {
////        //run the decider
////        FastqQCDecider decider = new FastqQCDecider();
////        decider.setProvenanceClient(extendedProvenanceClient);
////        decider.setWorkflowAccession(workflow.getSwAccession().toString());
////        run(decider, Arrays.asList(ArrayUtils.toArray(extraDeciderParams)));
//////        seqwareExecutor.deciderRunSchedule(FastqQCDecider.class.getCanonicalName(), workflow, extraDeciderParams);
////    }
//
//    private void run(OicrDecider decider, List<String> params) {
//        decider.setMetadata(metadata);
//        decider.setConfig(config);
//        decider.setParams(params);
//        decider.parse_parameters();
//        decider.init();
//        decider.do_test();
//        decider.do_run();
//        decider.clean_up();
//    }
//
//    @AfterClass
//    public void cleanUp() {
//        testEnv.shutdown();
//    }
//
//}
