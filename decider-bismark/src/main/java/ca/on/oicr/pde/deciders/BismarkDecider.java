package ca.on.oicr.pde.deciders;

import com.google.common.collect.Iterables;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author rtahir
 */
public class BismarkDecider extends OicrDecider {
    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;
    
    private final String [][] readMateFlags = {{"_R1_","1_sequence.txt",".1.fastq"},{"_R2_","2_sequence.txt",".2.fastq"}};    
    
    private String numOfThreads = "4";
    private String bismarkMemory   = "32";
    private String seedLength = "32";
    private String numMultiprocessingCores   = "2";
    private String maxMismatch = "1";
    private String queue = "";
    private String expectedOutputSam = "true";
    private String output_dir = "seqware-results"; 
    private String manualOutput = "false";
    private String templateType = "BS";
    private String output_prefix = "./";
    private String tmpDir = "tmp";
    
    private static final String FASTQ_GZ_METATYPE = "chemical/seq-na-fastq-gzip";
  

    public BismarkDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        this.setMetaType(Arrays.asList(FASTQ_GZ_METATYPE));
        this.setHeadersToGroupBy(Arrays.asList(FindAllTheFiles.Header.FILE_SWA));
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("manual-output", "Optional*. Set the manual output "
                + "either to true or false").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
                + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("tmp_dir", "Optional: Set the temp dir (Default: tmp)").withRequiredArg();
        //bismark
        parser.accepts("bismark_mem", "Optional: allocated memory required for bismark, default is 32.").withRequiredArg();
        //params
        parser.accepts("no_of_threads", "Optional: Bismark threads, default is 32.").withRequiredArg();
        parser.accepts("seed_length", "Optional: seed length").withRequiredArg();
        parser.accepts("max_mismatch_allowed", "Optional: maximum mismatch allowed").withRequiredArg();
        parser.accepts("no_of_multiprocessing_cores", "Optional: maximum mismatch allowed").withRequiredArg();
        // template type
        parser.accepts("temlate-type", "Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(FASTQ_GZ_METATYPE));
        this.setGroupingStrategy(FindAllTheFiles.Header.FILE_SWA);
        
        //bismark
        if (this.options.has("no_of_threads")) {
            this.numOfThreads = options.valueOf("no_of_threads").toString();
        }
        if (this.options.has("seed_length")) {
            this.seedLength = options.valueOf("seed_length").toString();
        }
        if (this.options.has("max_mismatch_allowed")) {
            this.maxMismatch = options.valueOf("max_mismatch_allowed").toString();
        }
        if (this.options.has("bismark_mem")) {
            this.bismarkMemory = options.valueOf("bismark_mem").toString();
        }

        ReturnValue rv = super.init();
        rv.setExitStatus(ReturnValue.SUCCESS);
        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            if (!options.hasArgument("template-type")) {
                Log.error("--template-type requires an argument, BS");
                rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                return rv;
            } else {
                this.templateType = options.valueOf("template-type").toString();
                if (!this.templateType.equals("BS")) {
                    Log.stderr("NOTE THAT ONLY BS template-type SUPPORTED, WE CANNOT GUARANTEE MEANINGFUL RESULTS WITH OTHER TEMPLATE TYPES");
                }
            }
        }

        if (this.options.has("manual-output")) {
            this.manualOutput = options.valueOf("manual_output").toString();
            Log.debug("Setting manual output, default is false and needs to be set only in special cases");
        }
        if (this.options.has("output-path")) {
            this.output_prefix = options.valueOf("output-path").toString();
            if (!this.output_prefix.endsWith("/")) {
                this.output_prefix += "/";
            }
        }

        if (this.options.has("output-folder")) {
            this.output_dir = options.valueOf("output-folder").toString();
        }

        return rv;
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (currentRV.getFiles().get(f).getMetaType().equals(FASTQ_GZ_METATYPE)) {
                        metatypeOK = true;
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }
            if (!metatypeOK) {
                continue; // Go to the next value
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(FindAllTheFiles.Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
        }

        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(FindAllTheFiles.Header.FILE_SWA.getTitle())).getGroupByAttribute();
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }
    
    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        // Filter the data of a different template type if filter is specified
        if (!this.templateType.equalsIgnoreCase(currentTtype)) {
            Log.warn("Excluding file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "] due to template type/geo_library_source_template_type = [" + currentTtype + "]");
            return false;
        }
        
        return super.checkFileDetails(returnValue, fm);
    }
    
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        if (filePaths.length != 2) {
            Log.error("This Decider supports only cases where we have only 2 files per lane, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }
        boolean haveFirstMate = false;
        boolean haveSecondMate = false;

        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                if (!haveFirstMate) {haveFirstMate = p.contains("R1");}
                if (!haveSecondMate){haveSecondMate= p.contains("R2");}
                    
                }
            }
        
        if (!haveFirstMate || !haveSecondMate) {
            Log.error("The Decider was not able to find both R1 and R2 fastq files for paired sequencing alignment, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }
        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }
    
    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    public ReturnValue customizeRun(WorkflowRun run) {
        ReturnValue rv = super.customizeRun(run);
        FileAttributes[] fas = run.getFiles();
//        int[] indexes = {0, 1};
        Set<String> sampleNames = new HashSet<>();

//        Set fqInputs_end1 = new HashSet();
//        Set fqInputs_end2 = new HashSet();
//        Set[] fqInputFiles = {fqInputs_end1, fqInputs_end2};
//        String fastq_inputs_end_1 = "";
//        String fastq_inputs_end_2 = "";
//        BeSmall currentBs = null;
        for (FileAttributes p : fas) {
            sampleNames.add(getRequiredAttribute(p, FindAllTheFiles.Header.SAMPLE_NAME));
        }
        
        String sampleName = null;
        if(sampleNames.size() != 1) {
            abortSchedulingOfCurrentWorkflowRun();
        } else {
            sampleName = Iterables.getOnlyElement(sampleNames);
        }
        
//            for (BeSmall bs : fileSwaToSmall.values()) {
//                if (!bs.getPath().equals(p)) {
//                    continue;
//                }
//
//                for (int i : indexes) {
//                    for (int j = 0; j < this.readMateFlags[i].length; j++) {
//                        if (p.toString().contains(this.readMateFlags[i][j])) {
//                            fqInputFiles[i].add(p);
//                            break;
//                        }
//                    }
//                }
//            currentBs = bs;
//            }
//        }
        // Refuse to continue if we don't have an object with metadta for one of the files
//        if (null == currentBs) {
//            Log.error("Was not able to retrieve fastq files for either one or two subsets of paired reads, not scheduling current workflow run");
//            this.abortSchedulingOfCurrentWorkflowRun();
//        }
        String fastq_inputs_end_1 = null;
        String fastq_inputs_end_2 = null;
        for (FileAttributes fa : fas) {
            if (fa.getPath().contains("R1")) {
                fastq_inputs_end_1 = fa.getPath();
            } else if (fa.getPath().contains("R2")) {
                fastq_inputs_end_2 = fa.getPath();
            } else {
                abortSchedulingOfCurrentWorkflowRun();
            }
        }

//        // Format input strings
//        if (fqInputFiles[0].size() == 0 || fqInputFiles[1].size() == 0) {
//            Log.error("Was not able to retrieve fastq files for either one or two subsets of paired reads, not scheduling current workflow run");
//            this.abortSchedulingOfCurrentWorkflowRun();
//        } else {
//            fastq_inputs_end_1 = _join(",", fqInputFiles[0]);
//            fastq_inputs_end_2 = _join(",", fqInputFiles[1]);
//        }

//        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        run.addProperty("r1_fastq_file", fastq_inputs_end_1);
        run.addProperty("r2_fastq_file", fastq_inputs_end_2);
        run.addProperty("max_mismatch_allowed", this.maxMismatch);
        run.addProperty("seed_length", this.seedLength);
        run.addProperty("no_of_threads", this.numOfThreads);
        run.addProperty("no_of_multiprocessing_cores", this.numMultiprocessingCores);
        run.addProperty("bismark_mem", this.bismarkMemory);
        run.addProperty("output_dir", this.output_dir);
        run.addProperty("template_type", this.templateType);
        run.addProperty("output_prefix", this.output_prefix);
        run.addProperty("manual_output", this.manualOutput);
        run.addProperty("tmp_dir", this.tmpDir);
        run.addProperty("output_filename_prefix", sampleName);
        
        if (!this.queue.isEmpty()) {
            run.addProperty("queue", this.queue);
        }
        

        return rv;
    }

   //Join function
   public static String _join(String separator, Set items) {
       StringBuilder result = new StringBuilder();
       Iterator myItems = items.iterator();
       while(myItems.hasNext()) {
          if (result.length() > 0)
              result.append(separator);

          result.append(myItems.next().toString());
       }

    return result.toString();
    }
   
    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(BismarkDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class BeSmall {

    private Date date = null;
    private String iusDetails = null;
    private String groupByAttribute = null;
    private String tissueType = null;
    private String path = null;
    private String tubeID = null;
    private String groupID = null;
    private String groupDescription = null;
    private String RGLB;
    private String RGPU;
    private String RGSM;
    //private String RGPL;
    private String ius_accession;
    private String sequencer_run_name;
    private String barcode;
    private String lane;

    public BeSmall(ReturnValue rv) {
        try {
            this.date = format.parse(rv.getAttribute(FindAllTheFiles.Header.PROCESSING_DATE.getTitle()));
        } catch (ParseException ex) {
            Log.error("Bad date!", ex);
            ex.printStackTrace();
        }

        FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
        this.tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
        this.tubeID = fa.getLimsValue(Lims.TUBE_ID);
        if (null == this.tubeID || this.tubeID.isEmpty()) {
            this.tubeID = "NA";
        }
        this.groupID = fa.getLimsValue(Lims.GROUP_ID);
        if (null == this.groupID || this.groupID.isEmpty()) {
            this.groupID = "NA";
        }
        this.groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
        if (null == this.groupDescription || this.groupDescription.isEmpty()) {
            this.groupDescription = "NA";
        }

        this.lane = fa.getLane().toString();
        this.RGLB = fa.getLibrarySample();
        this.RGPU = fa.getSequencerRun() + "_" + this.lane + "_" + fa.getBarcode();
        this.RGSM = fa.getDonor() + "_" + this.tissueType;
        if (!this.groupID.equals("NA")) {
            this.RGSM = this.RGSM + "_" + this.groupID;
        }

        this.iusDetails = this.RGLB + this.RGPU + rv.getAttribute(FindAllTheFiles.Header.FILE_SWA.getTitle());
        this.ius_accession = rv.getAttribute(FindAllTheFiles.Header.IUS_SWA.getTitle());
        this.sequencer_run_name = fa.getSequencerRun();
        this.barcode = fa.getBarcode();

        StringBuilder gba = new StringBuilder(fa.getDonor());
        gba.append(":").append(fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE));
        gba.append(":").append(this.ius_accession);

        String trs = fa.getLimsValue(Lims.TARGETED_RESEQUENCING);
        if (null != trs && !trs.isEmpty()) {
            gba.append(":").append(trs);
        }

        this.groupByAttribute = gba.toString();
        this.path = rv.getFiles().get(0).getFilePath() + "";
    }

    public Date getDate() {
        return this.date;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public String getGroupByAttribute() {
        return this.groupByAttribute;
    }

    public void setGroupByAttribute(String groupByAttribute) {
        this.groupByAttribute = groupByAttribute;
    }

    public String getTissueType() {
        return this.tissueType;
    }

    public String getIusDetails() {
        return this.iusDetails;
    }

    public void setIusDetails(String iusDetails) {
        this.iusDetails = iusDetails;
    }

    public String getPath() {
        return this.path;
    }

    public String getTubeId() {
        return this.tubeID;
    }

    public String getGroupID() {
        return this.groupID;
    }

    public String getGroupDescription() {
        return this.groupDescription;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getRGLB() {
        return RGLB;
    }

    public String getRGPU() {
        return RGPU;
    }

    public String getRGSM() {
        return RGSM;
    }

    public String getIus_accession() {
        return ius_accession;
    }

    public String getSequencer_run_name() {
        return sequencer_run_name;
    }

    public String getBarcode() {
        return barcode;
    }

    public String getLane() {
        return lane;
    }
}
}
