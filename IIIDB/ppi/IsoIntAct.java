package edu.nchu.syslab.ppi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.ResourceBundle;
import java.util.Set;

import edu.nchu.syslab.rnaseq.BuildExpMatrix;
import edu.nchu.syslab.rnaseq.ExpProfile;
import edu.usc.zhoulab.common.MathUtil;
import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.Enrich;
import edu.usc.zhoulab.common.annot.GeneAnnot;

public class IsoIntAct { 
	static ResourceBundle bundler = ResourceBundle.getBundle("common");
	static String DB_URL = bundler.getString("db.url");
	static String DB_USER = bundler.getString("db.user");
	static String DB_PASS = bundler.getString("db.passwd");
	
	//static String hg18_to_uniprot =	"/home/jimliu/ProcessUniprot/2011Mar/hg18_to_uniprot.txt";
	static String hg18_to_uniprot = null;
	static Connection m_conn = null;
	static String option = null;
	static String folderName = null;
	static boolean test_all_exp_flag = true;
	static String Tid_fileName = null;
	static String RatioMap_fileName = null;
	static String HomologyRatioMap_fileName = null;
	static String validation_option = null;
	static String GSPGSN_folderName = null;
	
	static String m_suffix = null;
	static boolean GSP_one2one_flag = true;
	static boolean TEST_one2one_flag = true;
	static double LR_th = 0.0;
	static String PredictionGO_fileName = null;
	static String PfamFile_fileName = null;
	static boolean m_debug_flag = true;
	static double sd_th = 0.1;
	static double corr_th = 0.0;
	
	//static Map<String, String> m_tidInt2pidInt = null;
	//static Set<String> m_pidInt = null;
	static Map<String, String> m_Known2geneidsMap = null;
	//static Map<String, Set<String>> m_Geneid2knownsMap = null;
	static Map<String, Set<String>> m_goid2tid = null;
	static Map<String, Set<String>> m_prot2goname = null;
	static Set<String> HomologyTid_set = null;
	Random random = new Random();
	
	static Set<String> Train_GSP_set = null;
	static Set<String> single_Train_GSP_set = null;
	static Set<String> Test_GSP_set = null;
	static Set<String> single_Test_GSP_set = null;
	static Set<String> Test_GSN_set = null;
	static Set<String> Train_GSN_set = null;
	static GSNset GSN_set = null;
	static List<String> tid_list = null;
	
	static boolean auto_getMaxMinLR_flag = true;
	static Double max_LR = 0.0;
	static Double min_LR = 0.0;
	
	static String m_pipeline = "Jim_pipeline";
	
	public static void main(String[] args) 
	throws Exception
	{
		if (args.length != 16) {
			System.out.println("Usage: IsoIntAct option" +
					" folderName" +
					" Tid_fileName" +
					" RatioMap_fileName m_suffix" +
					" GSP_one2one_flag test_all_exp_flag LR_th" +
					" PredictionGO_fileName" +
					" PfamFile_fileName" +
					" HomologyRatioMap_fileName" +
					" validation_option" +
					" auto_getMaxMinLR_flag max_LR min_LR" + 
					" GSPGSN_folderName");
			// outCoExp,
			// ENCODE_SRP000228_hg18_known_GM sel_tids.log buildRatioMap_n.log 0
			// true true 0.0
			// /home/jimliu/RNA-seq/isofrom_func/go2tid_p0.5_Oct2.txt 
			// /home/jimliu/RNA-seq/studyPfam/tid2pfam_f5.txt
			// /home/jimliu/RNA-seq/homology/knownInt_1e-6_0.1.txt
			// all CoExp DoEnrich SSBP Ortholog
			// /home/jimliu/ProcessUniprot/2011Mar/hg18_to_uniprot.txt;	
			// false 6316.0 0.0
			//manual set maximum
			//Double max_LR = 6316.0;
			//Double min_LR = 0.0;
			//SetA
		}
		option = args[0];
		folderName = args[1];
		Tid_fileName = args[2];
		RatioMap_fileName = args[3];
		m_suffix = args[4];
		GSP_one2one_flag = Boolean.parseBoolean(args[5]);
		test_all_exp_flag = Boolean.parseBoolean(args[6]);
		LR_th = Double.parseDouble(args[7]);
		PredictionGO_fileName = args[8];
		PfamFile_fileName = args[9];
		HomologyRatioMap_fileName = args[10];
		validation_option = args[11];
		auto_getMaxMinLR_flag = Boolean.parseBoolean(args[12]);
		max_LR = Double.parseDouble(args[13]); 
		min_LR = Double.parseDouble(args[14]);
		GSPGSN_folderName = args[15];
		System.out.println("option " + option);
		System.out.println("folderName " + folderName);
		System.out.println("Tid_fileName " + Tid_fileName);
		System.out.println("RatioMap_fileName " + RatioMap_fileName);
		System.out.println("m_suffix " + m_suffix);
		System.out.println("GSP_one2one_flag " + GSP_one2one_flag);
		System.out.println("test_all_exp_flag " + test_all_exp_flag);
		System.out.println("LR_th " + LR_th);
		System.out.println("PredictionGO_fileName " + PredictionGO_fileName);
		System.out.println("PfamFile_fileName " + PfamFile_fileName);
		System.out.println("HomologyRatioMap_fileName " + HomologyRatioMap_fileName);
		System.out.println("validation_option " + validation_option);
		System.out.println("hg18_to_uniprot " + hg18_to_uniprot);
		System.out.println("auto_getMaxMinLR_flag " + auto_getMaxMinLR_flag);
		System.out.println("max_LR " + max_LR);
		System.out.println("min_LR " + min_LR);
		System.out.println("GSPGSN_folderName " + GSPGSN_folderName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new IsoIntAct();
	}	
	
	public IsoIntAct() throws Exception {
		tid_list = FileUtil.getFileAsList(Tid_fileName);
		
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true);
		m_Known2geneidsMap = GeneAnnot.getKnown2geneid(m_conn);
		//m_Geneid2knownsMap = GeneAnnot.getGeneid2knowns(m_conn);
		
		//GSN_set = getGSNsetIntActIso("uniprot2go_0405_2011_BP", ""+GeneAnnot.TAX_HUMAN);
		//GSN_set = getGSNsetIntActIso("uniprot2go_0405_2011", ""+GeneAnnot.TAX_HUMAN);
		
		if (option.equals("outputGSN")){
			outputGSN();
		}
		
		/**Build 3 training/test sets:
		*  SetA. test set: GSP/GSN with group 1,  training set: GSP/GSN with group 2+3.
		*  SetB. test set: GSP/GSN with group 2,  training set: GSP/GSN with group 1+3.
		*  SetC. test set: GSP/GSN with group 3,  training set: GSP/GSN with group 1+2.
		*/
		
		if (option.equals("outputTrainTestSet")){
			Tid_List tid_list = new Tid_List();
			//tid_list = getTid_group_lsit(tid_list);
			tid_list.tid_list_group1 = FileUtil.getFileAsList("../tids_group1");
			tid_list.tid_list_group2 = FileUtil.getFileAsList("../tids_group2");
			tid_list.tid_list_group3 = FileUtil.getFileAsList("../tids_group3");
			getGSPGSN_group(tid_list);
			return;
		} else {
			
			// GSPGSN_folderName = /home/jimliu/IntAct/GSN_4_new2
			Train_GSP_set = FileUtil.getFileAsSet(GSPGSN_folderName	+ "/Train_GSP_set");
			single_Train_GSP_set = getSingle_tid_set(Train_GSP_set);
			Train_GSN_set = FileUtil.getFileAsSet(GSPGSN_folderName + "/Train_GSN_set");

			Test_GSP_set = FileUtil.getFileAsSet(GSPGSN_folderName + "/Test_GSP_set");
			single_Test_GSP_set = getSingle_tid_set(Test_GSP_set);
			Test_GSN_set = FileUtil.getFileAsSet(GSPGSN_folderName + "/Test_GSN_set");
		
		}
		
		/*if (option.equals("buildRatioMap")) {
			Set<Pair> GSP_set = getPairGSP(Train_GSP_set);
			buildRatioMap(GSP_set, RatioMap_fileName);
		}*/
		
		//TODO: option of calculate LR
		if (option.equals("outCoExp")) {
			ExpProfile transExpProfile = prepareExpProfile(tid_list);
			transExpProfile.outFile("transcript_exp.log");
			outCorr(transExpProfile, Train_GSP_set, Train_GSN_set);
		}
		if (option.equals("outDoEnrich")) {
			outRatioD(Train_GSP_set, Train_GSN_set);
		}
		if (option.equals("outSSBP")){
			outSSBP(Train_GSP_set, Train_GSN_set);
		}
		if (option.equals("outOrtholog")){
			outClass(Train_GSP_set, Train_GSN_set);
		}
		if (option.equals("outAllLR")){
			outRatioD(Train_GSP_set, Train_GSN_set);
			outSSBP(Train_GSP_set, Train_GSN_set);
			outClass(Train_GSP_set, Train_GSN_set);
		}
		
		//TODO: option of validation (PrecisionRecall)
		if (option.equals("validation")){
			List<ExpProfile> Profile_list = prepareExpProfile_list(tid_list);
			List<Map<String, Double>> CoExpLR_map_list = getLRmapList();
			Map<String, Double> DoEnrichLR_map = getLRmap("do_enrich_LR.log");
			Map<String, Double> SSBPLR_map = getLRmap("SSBP_LR.log");
			Map<String, Double> Ortholog_map = getLRmap("orthologs_LR.log");
			validation(Profile_list, CoExpLR_map_list, DoEnrichLR_map, SSBPLR_map, Ortholog_map);
			outPrecisionRecall(validation_option + "_interaction.log");
		}
		if (option.equals("PrecisionRecall")){
			outPrecisionRecall(validation_option + "_interaction.log");
		}

	}
	
	void outputGSN() throws Exception {
		FileUtil.outputCollectionAsFile(GSN_set.tids_nucleus, "nucleus_tids");
		FileUtil.outputCollectionAsFile(GSN_set.tids_membrane, "membrane_tids");
	}

	Tid_List getTid_group_lsit(Tid_List GSP_list) throws Exception {
		List<String> all_tid_list = FileUtil.getFileAsList(Tid_fileName);
		debug_println(" total  #isoforms: " + all_tid_list.size());
		
		GSP_list.tid_list_group1 = getRandomList(all_tid_list, all_tid_list.size()/3);
		debug_println(" group1 #isoforms: " + GSP_list.tid_list_group1.size());

		GSP_list.tid_list_group2 = getRandomList(all_tid_list, all_tid_list.size()/3, GSP_list.tid_list_group1);
		debug_println(" group2 #isoforms: " + GSP_list.tid_list_group2.size());
		
		GSP_list.tid_list_group3 = getTidList(all_tid_list, GSP_list.tid_list_group1, GSP_list.tid_list_group2);
		debug_println(" group3 #isoforms: " + GSP_list.tid_list_group3.size());
		return GSP_list;
	}
	
	List<String> getRandomList(List<String> tid_list, int size) {
		List<String> result = new ArrayList<String>();
		while (result.size() < size) {
			int i = random.nextInt(tid_list.size());
			String tid = tid_list.get(i);
			if (result.contains(tid))
				continue;
			result.add(tid);
		}
		return result;
	}
	
	List<String> getRandomList(List<String> tid_list, int size, List<String> tid_list1) {
		List<String> result = new ArrayList<String>();
		while (result.size() < size) {
			int i = random.nextInt(tid_list.size());
			String tid = tid_list.get(i);
			if (tid_list1.contains(tid))
				continue;
			if (result.contains(tid))
				continue;
			result.add(tid);
		}
		return result;
	}
	
	List<String> getTidList(List<String> tid_list, List<String> tid_list1, List<String> tid_list2) {
		List<String> result = new ArrayList<String>();
		for (String tid : tid_list) {
			if (tid_list1.contains(tid) || tid_list2.contains(tid))
				continue;
			if (result.contains(tid))
				continue;
			result.add(tid);
		}
		return result;
	}

	void getGSPGSN_group(Tid_List tid_list) throws Exception {
		Set<String> GSP_IntAct_set = getGSPsetIntActIso(m_conn, "/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		Set<String> single_GSP_set = getSingle_tid_set(GSP_IntAct_set);
		
		debug_println("output group1");
		outputGSPGSNfile("group1", tid_list.tid_list_group1,
				//intact_set.GSN_tid2tid_group1, intact_set.GSP_tid2tid_group1,
				GSP_IntAct_set, single_GSP_set);
		
		debug_println("output group2");
		outputGSPGSNfile("group2", tid_list.tid_list_group2,
				//intact_set.GSN_tid2tid_group2, intact_set.GSP_tid2tid_group2,
				GSP_IntAct_set, single_GSP_set);
		
		debug_println("output group3");
		outputGSPGSNfile("group3", tid_list.tid_list_group3,
				//intact_set.GSN_tid2tid_group3, intact_set.GSP_tid2tid_group3,
				GSP_IntAct_set, single_GSP_set);
	}
	
	void outputGSPGSNfile(String groupName, List<String> tid_group_list,
			//Set<String> GSN_tid2tid_set, Set<String> GSP_tid2tid_set,
			Set<String> GSP_IntAct_set, Set<String> single_GSP_set) throws Exception {
		Set<String> GSN_tid2tid_set = new HashSet<String>();
		Set<String> GSP_tid2tid_set = new HashSet<String>();
		int num_pos = 0;
		int num_neg = 0;
		for(int i=0; i < tid_group_list.size(); i++) {
			/*if (i % 1000 == 0)
				debug_println("#tid "+ i);*/
			String tid1 = tid_group_list.get(i);
			for(int j=0; j < tid_group_list.size(); j++) {
				String tid2 = tid_group_list.get(j);
				if (tid1.equals(tid2)){
					continue;
				}
				String tid2tid = tid1 + ":" + tid2;
				if (GSP_IntAct_set.contains(tid2tid)){
					GSP_tid2tid_set.add(tid2tid);
					num_pos++;
				}else if (single_GSP_set.contains(tid1) && single_GSP_set.contains(tid2)
							&& ((GSN_set.tids_membrane.contains(tid1) && GSN_set.tids_nucleus.contains(tid2)) || 
								(GSN_set.tids_membrane.contains(tid2) && GSN_set.tids_nucleus.contains(tid1)))) {
					num_neg++;
					GSN_tid2tid_set.add(tid2tid);
				}
			}
		}
		Set<String> single_GSP_tid_set = getSingle_tid_set(GSP_tid2tid_set);
		Set<String> single_GSN_tid_set = getSingle_tid_set(GSN_tid2tid_set);
		FileUtil.outputCollectionAsFile(tid_group_list, "tids_" + groupName);
		FileUtil.outputCollectionAsFile(single_GSP_tid_set, "GSP_isoforms_" + groupName);
		FileUtil.outputCollectionAsFile(GSP_tid2tid_set, "GSP_interact_" + groupName);
		FileUtil.outputCollectionAsFile(single_GSN_tid_set, "GSN_isoforms_" + groupName);
		FileUtil.outputCollectionAsFile(GSN_tid2tid_set, "GSN_interact_" + groupName);
		debug_println(" " + groupName + " cover #tids: " + single_GSP_tid_set.size());
		debug_println("    GSP cover #isoforms: " + single_GSP_tid_set.size());
		debug_println("    GSP #interactions: " + num_pos);
		debug_println("    GSN cover #isoforms: " + single_GSN_tid_set.size());
		debug_println("    GSN #interactions: " + num_neg);
		debug_println("    GSP/GSN: " + ((double)num_pos/num_neg));
	}

	static class Tid_List {
		List<String> tid_list_group1 = new ArrayList<String>();
		List<String> tid_list_group2 = new ArrayList<String>();
		List<String> tid_list_group3 = new ArrayList<String>();
	}
	
	ExpProfile prepareExpProfile(List<String> tid_list) throws Exception {
		new BuildExpMatrix(folderName, m_conn);
		List<String> sampleList = getSampleList(folderName);
		Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix.getTransMap(sampleList);
		ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(sampleList, sample2expMap);
		debug_println("transExpProfile: " + transExpProfile);
		transExpProfile.removeLowSD(sd_th); // sd_th 0.1			
		transExpProfile.applyGeneList(tid_list);
		return transExpProfile;
	}
	
	public static List<String> getSampleList(String folderName)
			throws IOException {
		List<String> SampleList = new ArrayList<String>();
		SampleList = FileUtil.getFilesInDir(folderName, "SRX.*");
		if (SampleList.size()==0){
			SampleList = FileUtil.getFilesInDir(folderName, "Myers.*");
		}
		if (SampleList.size()==0){
			SampleList = FileUtil.getFilesInDir(folderName, "SRR.*");
		}
		if (SampleList.size()==0){
			SampleList = FileUtil.getFilesInDir(folderName, "ERX.*");
		}
		System.out.println(SampleList);
		return SampleList;
	}
	
	static Set<String> getGSPsetIntActIso(Connection conn, String fileName, boolean one2one) throws Exception{
		debug_println("getGSPsetIntActIso()");
		debug_println("  hg18_to_uniprot " + hg18_to_uniprot);
		Set<String> known_set = new HashSet<String>();
		Set<String> knownInt = new HashSet<String>();
		Set<String> proteinInt = new HashSet<String>();
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		//m_tidInt2pidInt = new HashMap<String, String>();
		//m_pidInt = new HashSet<String>();
		
		List<String> lines = FileUtil.getFileAsList(fileName);
		Set<String> all_proteins_set = new HashSet<String>();	
		Set<String> sel_proteins_set = new HashSet<String>();
		for (String line: lines) {
			if (line.startsWith("#") || !line.contains("taxid:9606(Human)")) {
				continue;
			}
			String[] items = line.split("\t");
			String protein1 = items[0];
			String protein2 = items[1];
			if (!protein1.startsWith("uniprotkb:")) {
				//debug_println("  protein1: "+ protein1);
				continue;
			}
			if (!protein2.startsWith("uniprotkb:")) {
				//debug_println("  protein2: "+ protein2);
				continue;
			}
			{
				protein1 = protein1.replace('|', ' ');
				String[] items1 = protein1.split(" ");
				protein1 = items1[0].replace("uniprotkb:", "");
				protein2 = protein2.replace('|', ' ');
				String[] items2 = protein2.split(" ");
				protein2 = items2[0].replace("uniprotkb:", "");
			}
			proteinInt.add(protein1+":"+protein2);
			proteinInt.add(protein2+":"+protein1);
			
			all_proteins_set.add(protein1);
			all_proteins_set.add(protein2);
			
			if (protein1.equals(protein2)) {
				continue;
			}
			Set<String> knowns1 = protein2knowns.get(protein1);
			Set<String> knowns2 = protein2knowns.get(protein2);
			//debug_println("  protein1 "+ protein1 +" -> knowns "+ knowns1);
			if (knowns1 == null && protein1.contains("-")) {
				protein1 = protein1.substring(0, protein1.indexOf('-'));
				knowns1 = protein2knowns.get(protein1);
				//debug_println("    "+ protein1 +" -> knowns "+ knowns1);
			}
			if (knowns1 == null) {
				continue;
			}
			if (knowns2 == null && protein2.contains("-")) {
				protein2 = protein2.substring(0, protein2.indexOf('-'));
				knowns2 = protein2knowns.get(protein2);
			}
			if (knowns2 == null) {
				continue;
			}
			
			if (one2one) {
				if (knowns1.size() == 1 && knowns2.size() == 1) { // check one-to-one
					sel_proteins_set.add(protein1);
					sel_proteins_set.add(protein2);
					known_set.addAll(knowns1);
					known_set.addAll(knowns2);
					for (String g1 : knowns1) {
						for (String g2 : knowns2) {
							knownInt.add(g1 + ":" + g2);
							knownInt.add(g2 + ":" + g1);
							//m_tidInt2pidInt.put(g1 + ":" + g2, protein1 +":" + protein2);
							//m_tidInt2pidInt.put(g2 + ":" + g1, protein2 +":" + protein1);
							//m_pidInt.add(protein1 +":" + protein2);
							//m_pidInt.add(protein2 +":" + protein1);
						}
					}
				}
			} else {
				if (knowns1.size() >= 1 && knowns2.size() >= 1) {
					sel_proteins_set.add(protein1);
					sel_proteins_set.add(protein2);
					known_set.addAll(knowns1);
					known_set.addAll(knowns2);
					for (String g1 : knowns1) {
						for (String g2 : knowns2) {
							knownInt.add(g1 + ":" + g2);
							knownInt.add(g2 + ":" + g1);
							// TODO: build the map: tid-pair -> protein-pair
							//m_tidInt2pidInt.put(g1 + ":" + g2, protein1 +":" + protein2);
							//m_tidInt2pidInt.put(g2 + ":" + g1, protein2 +":" + protein1);
							//m_pidInt.add(protein1 +":" + protein2);
							//m_pidInt.add(protein2 +":" + protein1);
						}
					}
				}
			}
		}
		debug_println(fileName);
		debug_println("one2one "+ one2one);
		debug_println("proteinInt.size()/2 "+ (proteinInt.size()/2)); //total human PPIs: 
		debug_println("all_proteins_set "+ all_proteins_set.size()); //human proteins: 
		debug_println("sel_proteins_set "+ sel_proteins_set.size()); //human proteins with one-to-one mapping to knownGene: 
		debug_println("knownInt.size()/2 "+ (knownInt.size()/2)); //knownGene Interactions with one-to-one mapping to knownGene:
		debug_println("known_set "+ (known_set.size())); //knownGenes in IntAct database:
		return knownInt;
	}
	
	static class GSNset {
		Set<String> tids_nucleus = new HashSet<String>();
		Set<String> tids_both = new HashSet<String>();
		Set<String> tids_membrane = new HashSet<String>();
		//Set<String> geneids_nucleus = new HashSet<String>();
		//Set<String> geneids_membrane = new HashSet<String>();
	}
	
	/*GSNset getGSNsetIntActIso(Map<String, String> Known2geneids,
			Map<String, Set<String>> getGeneid2knowns) {
	debug_println("getGSNsetIntActIso()");
		Set<String> geneids = new HashSet<String>();
		for (String tid: Known2geneids.keySet()) {
			String geneid = Known2geneids.get(tid);
			geneids.add(geneid);
		}
		GSNset gsn = new GSNset();
		
		for (String g: geneids){
			Set<String> tids = getGeneid2knowns.get(g);
			Set<String> gos = Enrich.getGOs(g);
			if (gos==null)
				continue;
			if (gos.contains("GO:0005634 nucleus")){
				gsn.tids_nucleus.addAll(tids);
				gsn.tids_both.addAll(tids);
				gsn.geneids_nucleus.add(g);
			}
			if (gos.contains("GO:0005886 plasma membrane")){
				gsn.tids_membrane.addAll(tids);
				gsn.geneids_membrane.add(g);
			}
		}
		int num_both = 0;
		if (gsn.tids_both.retainAll(gsn.tids_membrane)){
			num_both = gsn.tids_both.size();
		}
		gsn.tids_membrane.removeAll(gsn.tids_both);
		gsn.tids_nucleus.removeAll(gsn.tids_both);
		int num_neg_iso = gsn.tids_nucleus.size() * gsn.tids_membrane.size();
		debug_println("# transcripts with nucleus: "+ gsn.tids_nucleus.size()); 
		debug_println("# transcripts with plasma membrane: "+ gsn.tids_membrane.size()); 
		debug_println("# transcripts with both nucleus and plasma membrane: "+ num_both); 
		debug_println("# negative interactions of isoforms: "+ num_neg_iso); 
		debug_println("# geneids_nucleus: "+ gsn.geneids_nucleus.size()); 
		debug_println("# geneids_membrane: "+ gsn.geneids_membrane.size()); 
		return gsn;
	}*/
	
	/**TODO: use database GO to build GSN
	 */
	GSNset getGSNsetIntActIso(String table, String tax) throws Exception {
		debug_println("getGSNsetIntActIso() table "+ table);	
		debug_println("  hg18_to_uniprot " + hg18_to_uniprot);
		debug_println("  tax "+ tax);
		
		GSNset gsn = new GSNset();
		
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		Map<String, Set<String>> tid2goname = new HashMap<String, Set<String>>();
		m_goid2tid = new HashMap<String, Set<String>>();

		String evi_condition = "and (evi = 'EXP' or evi = 'IDA' or evi = 'TAS' or evi = 'IC')";
		/*
		if (GSPGSN_folderName.equals("EXP_IDA")) {
			evi_condition = "and (evi = 'EXP' or evi = 'IDA')";
		}
		if (GSPGSN_folderName.equals("EXP_IC")) {
			evi_condition = "and (evi = 'EXP' or evi = 'IC' or evi = 'TAS' )";
		}
		if (GSPGSN_folderName.equals("IDA")) {
			evi_condition = "and (evi = 'IDA')";
		}
		if (GSPGSN_folderName.equals("TAS")) {
			evi_condition = "and (evi = 'TAS')";
		}
		if (GSPGSN_folderName.equals("IC_IDA")) {
			evi_condition = "and (evi = 'IC' or evi = 'IDA')";
		}*/
		debug_println("  evi_condition "+ evi_condition);
		PreparedStatement getTid = m_conn.prepareStatement("SELECT uniprot, go"
				+ "	FROM annot."+table
				+ " WHERE taxon = "+ tax +" "+evi_condition 
				+ " AND (go = 'GO:0005634' or go = 'GO:0005886')"
				+ " AND uniprot = ? ");
		for (String uniprot : protein2knowns.keySet()) {
			getTid.setString(1, uniprot);
			ResultSet rs = getTid.executeQuery();
			while (rs.next()) {
				String goid = rs.getString("go");
				Set<String> knowns = protein2knowns.get(uniprot);
				if (knowns == null) {
					continue;
				}
				for (String tid : knowns) {
					Set<String> gos = tid2goname.get(tid);
					if (gos == null) {
						gos = new HashSet<String>();
						tid2goname.put(tid, gos);
					}
					String goname = Enrich.getGoName(goid);
					if (goname != null) {
						tid2goname.get(tid).add(goname);
					}
				}
			}
		}
		
		for (String tid : tid_list){
			Set<String> gos = tid2goname.get(tid);
			if (gos == null) {
				continue;
			}
			if (gos.contains("GO:0005634 nucleus")){
				gsn.tids_nucleus.add(tid);
				gsn.tids_both.add(tid);
			}
			if (gos.contains("GO:0005886 plasma membrane")){
				gsn.tids_membrane.add(tid);
			}
		}
		int num_both = 0;
		if (gsn.tids_both.retainAll(gsn.tids_membrane)){
			num_both = gsn.tids_both.size();
		}
		gsn.tids_membrane.removeAll(gsn.tids_both);
		gsn.tids_nucleus.removeAll(gsn.tids_both);
		//int num_neg_iso = gsn.tids_nucleus.size() * gsn.tids_membrane.size();
		debug_println("# transcripts with nucleus: "+ gsn.tids_nucleus.size()); 
		debug_println("# transcripts with plasma membrane: "+ gsn.tids_membrane.size()); 
		//debug_println("# transcripts with both nucleus and plasma membrane: "+ num_both); 
		//debug_println("# negative interactions of isoforms: "+ num_neg_iso);
		
		String membrane_fileName = "/home/ting/IntAct/HttpUnit/result/HummPLoc2_result_membrane";
		String nucleus_fileName = "/home/ting/IntAct/HttpUnit/result/HummPLoc2_result_nucleus";
		Map<String, String> PLoc_membrane_map =  FileUtil.getFileAsMap(membrane_fileName);
		Map<String, String> PLoc_nucleus_map =  FileUtil.getFileAsMap(nucleus_fileName);
		for (String tid: PLoc_membrane_map.keySet()){
			if (PLoc_membrane_map.get(tid).contains("Nucleus") ||
					!PLoc_membrane_map.get(tid).contains("Cell membrane") ){
				gsn.tids_membrane.remove(tid);
			}
		}
		for (String tid: PLoc_nucleus_map.keySet()){
			if (PLoc_nucleus_map.get(tid).contains("Cell membrane") ||
					!PLoc_nucleus_map.get(tid).contains("Nucleus") ){
				gsn.tids_nucleus.remove(tid);
			}
		}
		int num_neg_iso = gsn.tids_nucleus.size() * gsn.tids_membrane.size();
		debug_println("# after remove not consistent result from PLoc");
		debug_println("# transcripts with nucleus: "+ gsn.tids_nucleus.size()); 
		debug_println("# transcripts with plasma membrane: "+ gsn.tids_membrane.size()); 
		debug_println("# transcripts with both nucleus and plasma membrane: "+ num_both); 
		debug_println("# negative interactions of isoforms: "+ num_neg_iso);
		return gsn;
	}
	
	
	Map<String, Set<String>> getTid2gonameUniProt(String table)	throws Exception {
		return getTid2gonameUniProt(table, ""+GeneAnnot.TAX_HUMAN);
	}
	
	Map<String, Set<String>> getTid2gonameUniProt(String table, String tax) 
	throws Exception {
		debug_println("getTid2gonameUniProt() table "+ table);
		debug_println("  hg18_to_uniprot " + hg18_to_uniprot);
		//table = uniprot2go_0405_2011
		debug_println("  tax "+ tax);

		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		Map<String, Set<String>> tid2goname = new HashMap<String, Set<String>>();
		m_goid2tid = new HashMap<String, Set<String>>();
		//m_prot2goname = new HashMap<String, Set<String>>();
	    //EXP, IDA, TAS, IC: the top level annotation.
	    //IPI, IMP, IGI, IEP: the second level annotation.

		//String evi_condition = "and (evi = 'TAS' or evi = 'IC')";
		//String evi_condition = "and (evi = 'IC')";
		
		//Select top level evidence code: EXP, IDA, TAS, IC 
		String evi_condition = "and (evi = 'EXP' or evi = 'IDA' or evi = 'TAS' or evi = 'IC')";
		
		//String evi_condition = " and (evi != 'IEA')";
		//String evi_condition = " and (evi = 'IDA')";
		//String evi_condition = " ";
		debug_println("  evi_condition "+ evi_condition);
		PreparedStatement getTid = m_conn.prepareStatement("SELECT uniprot, go"
				+ "	FROM annot."+table
				+ " WHERE taxon = "+ tax +" "+evi_condition 
				+ " AND uniprot = ? ");
		for (String uniprot : protein2knowns.keySet()) {
			getTid.setString(1, uniprot);
			ResultSet rs = getTid.executeQuery();
			while (rs.next()) {
				String goid = rs.getString("go");
				Set<String> knowns = protein2knowns.get(uniprot);
				if (knowns == null) {
					continue;
				}
				for (String tid : knowns) {
					Set<String> gos = tid2goname.get(tid);
					if (gos == null) {
						gos = new HashSet<String>();
						tid2goname.put(tid, gos);
					}
					String goname = Enrich.getGoName(goid);
					if (goname != null) {
						gos.add(goname);
					}
					Set<String> tids = m_goid2tid.get(goid);
					if (tids == null) {
						tids = new HashSet<String>();
						m_goid2tid.put(goid, tids);
					}
					tids.add(tid);
				}
			}
		}	
		debug_println("  tid2goname "+ tid2goname.size());
		debug_println("  m_goid2tid "+ m_goid2tid.size());
		return tid2goname;
	}
	
	public static Set<String> getSingle_tid_set(Set<String> tid2tid_set) throws Exception{
		Set<String> single_tid_set = new HashSet<String>();	
		for (String t: tid2tid_set){
			String[] items = t.split(":");
			single_tid_set.add(items[0]);
			single_tid_set.add(items[1]);
		}
		debug_println("single_tid_set size "+ single_tid_set.size());
		return single_tid_set;
	}
	
	//TODO: calculate LR
	void outCorr(ExpProfile profile, Set<String> Train_GSP_set, Set<String> Train_GSN_set) throws Exception {
		//Map<String, Double> Corr_map = new HashMap<String, Double>();
		int[] gsp = new int[21];
		int[] gsn = new int[21];
		int[] total = new int[21];
		debug_println("sd_th "+ sd_th);		
		List<String> idList = profile.idList;
		double[][] exp = profile.exp;
		for(int i=0; i < exp.length; i++) {
			if (i % 1000 == 0)
				debug_println("#tid "+ i);
			for(int j=i+1; j < exp.length; j++) {
				String tid1 = idList.get(i);
				String tid2 = idList.get(j);
				if (tid1.equals(tid2)){
					continue;
				}
				double corr = 0.0;
				if ((exp[i] != null)&&(exp[j] != null)) {
					corr = MathUtil.getPearson(exp[i], exp[j]);
				}
				int index = (int)(corr * 10) + 10;
				total[index]++;
				
				String PPI = tid1 + ":" + tid2;
				if (Train_GSP_set.contains(PPI)){
					gsp[index]++;					
				}
				if (Train_GSN_set.contains(PPI)){
					gsn[index]++;		
				}
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("co_exp_LR_"+m_suffix+".log"));
		int sum_gsp = MathUtil.getSum(gsp);
		int sum_gsn = MathUtil.getSum(gsn);
		int sum_total = MathUtil.getSum(total);
		bw.write("Index\tR\tGSP\tGSN\tTOTAL\tPr(R|GSP)\tPr(R|GSN)\tLR\t-");
		bw.write("\n");
		for (int i = 0; i < total.length; i++){
			double R = (double) (i-10) / 10;
			double p_gsp = (double) gsp[i] / sum_gsp;
			double p_gsn = (double) gsn[i] / sum_gsn;
			double LR = 0.0;
			if (p_gsn != 0){
				LR = p_gsp / p_gsn;
			}
			//Corr_map.put(Double.toString(R), LR);
			bw.write(i+"\t"+R+"\t"+gsp[i]+"\t"+gsn[i]+"\t"+total[i]+"\t"+p_gsp+"\t"+p_gsn+"\t"+LR);
			bw.write("\n");
		}	
		bw.write("-\tSUM\t"+sum_gsp+"\t"+sum_gsn+"\t"+sum_total);
		bw.close();
		//return Corr_map;
	}
	
	//D = Pr (di : dj | GSP) / ( Pr(di | GSP) * Pr(dj | GSP) )    di : dj >= 3

	static Map<String, Double> getRatioMap(String fileName) throws Exception{
		Map<String, Double> result = new HashMap<String, Double>();
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		while ((line = br.readLine()) != null) {
			String[] items = line.split("\t");
			if (items.length == 2) {
				result.put(items[0], Double.parseDouble(items[1]));
			}
		}
		br.close();
		debug_println("getRatioMap() RatioMap "+ result.size());
		return result;		
	}	
	
	static Map<String, Double> getHomologyRatioMap(String fileName) throws Exception{
		HomologyTid_set = new HashSet<String>();
		Map<String, Double> result = new HashMap<String, Double>();
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		double maxRatio = 0.0;
		while ((line = br.readLine()) != null) {
			String[] items = line.split("\t");
			if (items.length != 3)
				continue;
			HomologyTid_set.add(items[0]);
			HomologyTid_set.add(items[1]);
			if (result.get(items[0] + ":" + items[1]) == null) {
				result.put(items[0] + ":" + items[1], Double.parseDouble(items[2]));
			} else {
				maxRatio = result.get(items[0] + ":" + items[1]);
				if ( Double.parseDouble(items[2]) > maxRatio)
					maxRatio = Double.parseDouble(items[2]);
				result.put(items[0] + ":" + items[1], maxRatio);
			}
		}
		br.close();
		debug_println("getHomologyRatioMap() HomologyTid_set "+ HomologyTid_set.size());
		debug_println("getHomologyRatioMap() HomologyRatioMap "+ result.size());
		return result;		
	}
	
	void outRatioD(Set<String> Train_GSP_set, Set<String> Train_GSN_set) throws Exception {
		//Map<String, Double> RatioD_map = new HashMap<String, Double>();
		int[] gsp = new int[5];
		int[] gsn = new int[5];
		int[] total = new int[5];
		Map<String, Set<String>> tid2Pfam = new HashMap<String, Set<String>>();
		if (!PfamFile_fileName.equals("null")){
			tid2Pfam = getTid2PfamFile(PfamFile_fileName);
		}else{
			tid2Pfam = getTid2Pfam();
		}
		Set<String> tid_set = tid2Pfam.keySet();
		Map<String, Double> RatioMap = getRatioMap(RatioMap_fileName);
		int count = 0;
		for(String tid1: tid_set) {
			count++;
			/*if (count >= 10)
				break;*/
			if (count % 1000 == 0)
				debug_println("tid_set "+ count);
			for(String tid2: tid_set) {
				if (tid1.equals(tid2)){
					continue;
				}
				Set<String> domains_i = tid2Pfam.get(tid1);
				Set<String> domains_j = tid2Pfam.get(tid2);
				double maxRatioD = 0.0;
				double RatioD = 0.0;
				//double[] RatioD_array = new double[domains_i.size()*domains_j.size()+1];
				//int do_count = 0;
				for (String di: domains_i) {
					for (String dj: domains_j){
						/*
						if (di.equals(dj))
							continue; */
						//do_count++;
						String domain_pair = di + ":" + dj;
						if (RatioMap.get(domain_pair) != null)
							RatioD = RatioMap.get(domain_pair);
						if (maxRatioD < RatioD) {
							maxRatioD = RatioD;
						}
						//RatioD_array[do_count] = RatioD;
					}
				}
				//maxRatioD = MathUtil.getMax(RatioD_array);
				int index = 0;
				if (maxRatioD >= 10) index = 1;
				else if (maxRatioD >= 3) index = 2;
				else if (maxRatioD >= 2) index = 3;
				else if (maxRatioD >= 1) index = 4;
				total[index]++;
				String PPI = tid1 + ":" + tid2;
				if (Train_GSP_set.contains(PPI)){
					gsp[index]++;					
				}
				if (Train_GSN_set.contains(PPI)){
					gsn[index]++;		
				}	
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("do_enrich_LR.log"));
		int sum_gsp = MathUtil.getSum(gsp);
		int sum_gsn = MathUtil.getSum(gsn);
		int sum_total = MathUtil.getSum(total);
		bw.write("Index\tD\tGSP\tGSN\tTOTAL\tPr(D|GSP)\tPr(D|GSN)\tLR\t-");
		bw.write("\n");
		for (int i = 1; i < total.length; i++){
			int D = 0;
			if (i == 1)	D = 10;
			if (i == 2)	D = 3;
			if (i == 3)	D = 2;
			if (i == 4)	D = 1;
			double p_gsp = (double) gsp[i] / sum_gsp;
			double p_gsn = (double) gsn[i] / sum_gsn;
			double LR = 0.0;
			if (p_gsn != 0){
				LR = p_gsp / p_gsn;
			}
			//RatioD_map.put(Integer.toString(D), LR);
			bw.write(i+"\t"+D+"\t"+gsp[i]+"\t"+gsn[i]+"\t"+total[i]+"\t"+p_gsp+"\t"+p_gsn+"\t"+LR);
			bw.write("\n");
		}	
		bw.write("-\tpossible\t"+sum_gsp+"\t"+sum_gsn+"\t"+sum_total);
		bw.close();
		//return RatioMap;
	}
	
	void outSSBP(Set<String> Train_GSP_set, Set<String> Train_GSN_set) throws Exception {
		int[] gsp = new int[7];
		int[] gsn = new int[7];
		int[] total = new int[7];
		Map<String, Set<String>> Tid2goname = new HashMap<String, Set<String>>();
		if (!PredictionGO_fileName.equals("null")){
			Tid2goname = getTid2gonamePrediction(PredictionGO_fileName);
		}else{
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0723_2009");
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0728_2010");
			Tid2goname = getTid2gonameUniProt("uniprot2go_0405_2011_BP");
		}
		Set<String> tid_set = Tid2goname.keySet();
		int count = 0;
		for(String tid1: tid_set) {
			count++;
			/*if (count >= 10)
				break;*/
			if (count % 1000 == 0)
				debug_println("tid_set "+ count);
			if (Tid2goname.get(tid1) == null){
				continue;
			}
			for(String tid2: tid_set) {
				if (Tid2goname.get(tid2) == null){
					continue;
				}
				if (tid1.equals(tid2)){
					continue;
				}
				Set<String> gos_intersection = new HashSet<String>(Tid2goname.get(tid1));
				//Set<String> gos_i = Tid2goname.get(tid1);
				Set<String> gos_j = Tid2goname.get(tid2);
				gos_intersection.retainAll(gos_j);
				if (gos_intersection.size() == 0)
					continue;
				//System.out.println(tid1 + ":" + tid2);
				//System.out.println("gos_intersection.size() " + gos_intersection.size());
				int SSBP = 0;
				int[] SSBP_array = new int[gos_intersection.size()+1];
				int go_count = 0;
				for (String go: gos_intersection) {
					String[] items = go.split(" ");
					go = items[0];
					//System.out.print("go " + go + " ");
					if (m_goid2tid.get(go).size() > 0){
						go_count++;
						SSBP_array[go_count] = m_goid2tid.get(go).size();
						//System.out.println(m_goid2tid.get(go).size());
					}
				}
				SSBP = getMin(SSBP_array);
				int index = 0;
				if (SSBP < 10) index = 1;
				else if (SSBP < 50) index = 2;
				else if (SSBP < 100) index = 3;
				else if (SSBP < 500) index = 4;
				else if (SSBP < 1000) index = 5;
				else if (SSBP < 1000000) index = 6;
				//System.out.print(index + ":" + SSBP + "\n");
				total[index]++;
				String PPI = tid1 + ":" + tid2;
				if (Train_GSP_set.contains(PPI)){
					gsp[index]++;					
				}
				if (Train_GSN_set.contains(PPI)){
					gsn[index]++;		
				}
				/*if ((GSN_set.tids_membrane.contains(tid1) && 
					 GSN_set.tids_nucleus.contains(tid2)) ||
					(GSN_set.tids_membrane.contains(tid2) && 
					 GSN_set.tids_nucleus.contains(tid1))
					 ) {
					gsn[index]++;
				}*/
				
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("SSBP_LR.log"));
		int sum_gsp = MathUtil.getSum(gsp);
		int sum_gsn = MathUtil.getSum(gsn);
		int sum_total = MathUtil.getSum(total);
		bw.write("Index\tSSBP\tGSP\tGSN\tTOTAL\tPr(D|GSP)\tPr(D|GSN)\tLR\t-");
		bw.write("\n");
		for (int i = 1; i < total.length; i++){
			int SSBP = 0;
			if (i == 1)	SSBP = 10;
			if (i == 2)	SSBP = 50;
			if (i == 3)	SSBP = 100;
			if (i == 4)	SSBP = 500;
			if (i == 5)	SSBP = 1000;
			if (i == 6)	SSBP = 1000000;
			double p_gsp = (double) gsp[i] / sum_gsp;
			double p_gsn = (double) gsn[i] / sum_gsn;
			double LR = 0.0;
			if (p_gsn != 0){
				LR = p_gsp / p_gsn;
			}
			//RatioD_map.put(Integer.toString(D), LR);
			bw.write(i+"\t"+SSBP+"\t"+gsp[i]+"\t"+gsn[i]+"\t"+total[i]+"\t"+p_gsp+"\t"+p_gsn+"\t"+LR);
			bw.write("\n");
		}	
		bw.write("-\tpossible\t"+sum_gsp+"\t"+sum_gsn+"\t"+sum_total);
		bw.close();
		
	}	
	
	public static int getMin(int[] x) {
		int min = 10000000;
		for (int i = 0; i < x.length; i++) {
			if (x[i] == 0) 
				continue;
			if (x[i] < min) {
				min = x[i];
			}
		}
		return min;
	}

	
	/*static Map<String, Set<String>> getTid2gonameUniProt(String table) throws Exception {
		debug_println("getTid2gonameUniProt() table "+ table);
		
		m_goid2tid = new HashMap<String, Set<String>>();
		Map<String, Set<String>> tid2goname = new HashMap<String, Set<String>>();
		PreparedStatement getTid = m_conn.prepareStatement("SELECT k.name, k.proteinID, u.go "
				+ "	FROM ucsc_hg18.knownGene k, annot."+table+" u " 
				+ " WHERE k.proteinID = u.uniprot ");
		PreparedStatement getTid = m_conn.prepareStatement("SELECT k.name, k.proteinID, u.go"
				+ " FROM (SELECT * FROM annot."+table+" u" 
				+ " WHERE evi = 'EXP' or evi = 'IDA' or evi = 'TAS' or evi = 'IC') u, ucsc_hg18.knownGene k"
				+ " WHERE k.proteinID = u.uniprot");
		ResultSet rs = getTid.executeQuery();			
		while (rs.next()) {
			String tid = rs.getString("k.name");
			//String uniprot = rs.getString("k.proteinID");			
			String goid = rs.getString("u.go");
			if (Enrich.getGoName(goid) == null) {
				continue;
			}
			{
				Set<String> go = tid2goname.get(tid);
				if (go == null) {
					go = new HashSet<String>();					
					tid2goname.put(tid, go);
				}
				String goname = Enrich.getGoName(goid);
				if (goname != null) {
					go.add(goname);
				}
			}
			{
				Set<String> tids = m_goid2tid.get(goid);
				if (tids == null) {
					tids = new HashSet<String>();
					m_goid2tid.put(goid, tids);
				}
				tids.add(tid);
			}			
		}
		debug_println("getTid2gonameUniProt() tid2goname "+ tid2goname.size());
		debug_println("  m_goid2tid "+ m_goid2tid.size());
		return tid2goname;
	}*/
	
	void outClass(Set<String> Train_GSP_set, Set<String> Train_GSN_set) throws Exception {
		//Map<String, Double> RatioD_map = new HashMap<String, Double>();
		int[] gsp = new int[6];
		int[] gsn = new int[6];
		int[] total = new int[6];
		Map<String, Double> RatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);
		int count = 0;
		for(String tid1: HomologyTid_set) {
			count++;
			if (count % 1000 == 0)
				debug_println("tid_set "+ count);
			for(String tid2: HomologyTid_set) {
				if (tid1.equals(tid2)){
					continue;
				}
				double Class = 0.0;
				String PPI = tid1 + ":" + tid2;
				if (RatioMap.get(PPI) != null)
					Class = RatioMap.get(PPI);
				int index = 0;
				if (Class > 0.8) index = 5;
				else if (Class > 0.6) index = 4;
				else if (Class > 0.4) index = 3;
				else if (Class > 0.2) index = 2;
				else if (Class > 0.0) index = 1;
				total[index]++;

				if (Train_GSP_set.contains(PPI)){
					gsp[index]++;					
				}
				if (Train_GSN_set.contains(PPI)){
					gsn[index]++;		
				}
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("orthologs_LR.log"));
		int sum_gsp = MathUtil.getSum(gsp);
		int sum_gsn = MathUtil.getSum(gsn);
		int sum_total = MathUtil.getSum(total);
		bw.write("Index\tCLASS\tGSP\tGSN\tTOTAL\tPr(CL|GSP)\tPr(CL|GSN)\tLR\t-");
		bw.write("\n");
		for (int i = 1; i < total.length; i++){
			//double CLASS = (double) i / 10;
			double CLASS = 0.0;
			if (i == 1)	CLASS = 0.0;
			if (i == 2)	CLASS = 0.2;
			if (i == 3)	CLASS = 0.4;
			if (i == 4)	CLASS = 0.6;
			if (i == 5)	CLASS = 0.8;
			double p_gsp = (double) gsp[i] / sum_gsp;
			double p_gsn = (double) gsn[i] / sum_gsn;
			double LR = 0.0;
			if (p_gsn != 0){
				LR = p_gsp / p_gsn;
			}
			//RatioD_map.put(Integer.toString(D), LR);
			bw.write(i+"\t"+CLASS+"\t"+gsp[i]+"\t"+gsn[i]+"\t"+total[i]+"\t"+p_gsp+"\t"+p_gsn+"\t"+LR);
			bw.write("\n");
		}	
		bw.write("-\tpossible\t"+sum_gsp+"\t"+sum_gsn+"\t"+sum_total);
		bw.close();
		//return RatioMap;
	}
	
	
	Map<String, Double> getLRmap(String fileName) throws Exception {
		Map<String, Double> LR_map = new HashMap<String, Double>();
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		while ((line = br.readLine()) != null) {
			String[] items = line.split("\t");
			if (items.length == 8) {
				String key = items[1];
				Double LR = Double.parseDouble(items[7]);
				LR_map.put(key, LR);
			}
		}
		br.close();
		return LR_map;
	}
	
	//TODO: validation of all data
	List<ExpProfile> prepareExpProfile_list(List<String> tid_list) throws Exception {		
		List<ExpProfile> Profile_list = new ArrayList<ExpProfile>();
		if (m_pipeline != null) {
			List<String> datasets = FileUtil.getFileAsList("/home/jimliu/RNA-seq/Results/" + m_pipeline	+ "/dataset.lst");
			for (String d : datasets) {
				folderName = "/home/jimliu/RNA-seq/Results/" + m_pipeline + "/"	+ d;
				debug_println("RNA-seq " + folderName);
				new BuildExpMatrix(folderName, m_conn);
				List<String> sampleList = getSampleList(folderName);
				Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix.getTransMap(sampleList);
				ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(sampleList, sample2expMap);
				debug_println("transExpProfile " + d + ": " + transExpProfile);
				transExpProfile.removeLowSD(sd_th); // sd_th 0.1
				transExpProfile.applyGeneList(tid_list);
				Profile_list.add(transExpProfile);
			}
		} else {
			for (int i = 1; i <= 11; i++) {
				if (i == 6 || i == 8)
					continue;
				folderName = "/home/jimliu/RNA-seq/Results/d" + i + "_hg18_known";
				debug_println("RNA-seq " + folderName);
				new BuildExpMatrix(folderName, m_conn);
				List<String> sampleList = getSampleList(folderName);
				Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix.getTransMap(sampleList);
				ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(sampleList, sample2expMap);
				debug_println("transExpProfile d" + i + ": " + transExpProfile);
				transExpProfile.removeLowSD(sd_th); // sd_th 0.1
				transExpProfile.applyGeneList(tid_list);
				Profile_list.add(transExpProfile);
			}
		}
		debug_println("prepareExpProfile_list() " + Profile_list.size());
		return Profile_list;
	}
	
	List<Map<String, Double>> getLRmapList() throws Exception {
		List<Map<String, Double>> CoExpLR_map_list = new ArrayList<Map<String, Double>>();
		if (m_pipeline != null) {
			List<String> datasets = FileUtil.getFileAsList("/home/jimliu/RNA-seq/Results/" + m_pipeline	+ "/dataset.lst");
			for (String d : datasets) {
				folderName = "/home/jimliu/RNA-seq/Results/" + m_pipeline + "/"	+ d;
				debug_println("RNA-seq " + folderName);
				String fileName = "/home/ting/IntAct/CoExpLR_sd01/co_exp_LR_" +d+ ".log";
				Map<String, Double> LR_map = new HashMap<String, Double>();
				BufferedReader br = new BufferedReader(new FileReader(fileName));
				String line;
				while ((line = br.readLine()) != null) {
					String[] items = line.split("\t");
					if (items.length == 8) {
						String key = items[1];
						Double LR = Double.parseDouble(items[7]);
						LR_map.put(key, LR);
					}
				}
				br.close();
				debug_println("co_exp_LR_" + d);
				CoExpLR_map_list.add(LR_map);
			}
		} else {
		for (int i = 1; i <= 11; i++){
			if (i==6 || i==8) continue;
				String fileName = "/home/ting/IntAct/CoExpLR_sd01/co_exp_LR_d"+i+".log";
				Map<String, Double> LR_map = new HashMap<String, Double>();
				BufferedReader br = new BufferedReader(new FileReader(fileName));
				String line;
				while ((line = br.readLine()) != null) {
					String[] items = line.split("\t");
					if (items.length == 8) {
						String key = items[1];
						Double LR = Double.parseDouble(items[7]);
						LR_map.put(key, LR);
					}
				}
				br.close();
				debug_println("co_exp_LR_d"+i);
				CoExpLR_map_list.add(LR_map);
			}
		}
		debug_println("getLRmapList() " + CoExpLR_map_list.size());
		return CoExpLR_map_list;
	}
	
	double getExpMaxLR(int i, int j, List<ExpProfile> profile_list,
			List<Map<String, Double>> CoExpLR_map_list)	throws Exception {
		//double R = 0;
		double maxCoExp_LR = 0;
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			double R = 0;
			if (exp[i] != null && exp[j] != null)
				R = MathUtil.getPearson(exp[i], exp[j]);

			String R_str = Double.toString(R);
			if (R_str.startsWith("-")) {
				R_str = R_str.substring(0, 4);
			} else {
				R_str = R_str.substring(0, 3);
			}
			if (R_str.startsWith("-0.0")) {
				R_str = "0.0";
			}
			double CoExp_LR = 1.0;
			if (CoExpLR_map_list.get(net).get(R_str) != null) {
				CoExp_LR = CoExpLR_map_list.get(net).get(R_str);
			}
			if (maxCoExp_LR < CoExp_LR) {
				maxCoExp_LR = CoExp_LR;
			}
		}
		if (maxCoExp_LR == 0.0)
			maxCoExp_LR = 1.0;
		return maxCoExp_LR;
	}
	
	double getExpLR(int i, int j, List<ExpProfile> profile_list,
			List<Map<String, Double>> CoExpLR_map_list)	throws Exception {
		int net = 0;//net = 0 -> co_exp_LR_d2; net = 1 -> co_exp_LR_d3
		double[][] exp = profile_list.get(net).exp;
		double R = 0;
		if (exp[i] != null && exp[j] != null)
			R = MathUtil.getPearson(exp[i], exp[j]);

		String R_str = Double.toString(R);
		if (R_str.startsWith("-")) {
			R_str = R_str.substring(0, 4);
		} else {
			R_str = R_str.substring(0, 3);
		}
		if (R_str.startsWith("-0.0")) {
			R_str = "0.0";
		}
		double CoExp_LR = 1.0;
		if (CoExpLR_map_list.get(net).get(R_str) != null) {
			CoExp_LR = CoExpLR_map_list.get(net).get(R_str);
		}
		if (CoExp_LR == 0.0)
			CoExp_LR = 1.0;
		return CoExp_LR;
	}
	
	double getDomainLR(String tid1, String tid2, Map<String, Set<String>> tid2Pfam, 
			Map<String, Double> RatioMap, Map<String, Double> DoEnrichLR_map) throws Exception {
		double maxRatioD = 0.0;
		double RatioD = 0.0;
		if (tid2Pfam.get(tid1) != null && tid2Pfam.get(tid2) != null){
			Set<String> domains_i = tid2Pfam.get(tid1);
			Set<String> domains_j = tid2Pfam.get(tid2);
			//double[] RatioD_array = new double[domains_i.size()*domains_j.size()+1];
			int do_count = 0;
			for (String di: domains_i) {
				for (String dj: domains_j){
					if (di.equals(dj))
						continue;
					do_count++;
					String domain_pair = di + ":" + dj;
					if (RatioMap.get(domain_pair) != null)
						RatioD = RatioMap.get(domain_pair);
					if (maxRatioD < RatioD)
						maxRatioD = RatioD;
					//RatioD_array[do_count] = RatioD;
				}
			}
			//RatioD = MathUtil.getMax(RatioD_array);
		}
		String RatioD_str = "0";
		/*if (maxRatioD > 10) RatioD_str = "10";
		else if (maxRatioD > 3) RatioD_str = "3";
		else if (maxRatioD > 2) RatioD_str = "2";
		else if (maxRatioD > 1) RatioD_str = "1";*/
		if (maxRatioD >= 10) RatioD_str = "10";
		else if (maxRatioD >= 3) RatioD_str = "3";
		else if (maxRatioD >= 2) RatioD_str = "2";
		else if (maxRatioD >= 1) RatioD_str = "1";
		double DoEnrich_LR = 1.0;
		if (DoEnrichLR_map.get(RatioD_str)!=null){
			DoEnrich_LR = DoEnrichLR_map.get(RatioD_str);
		}
		return DoEnrich_LR;
	}

	double getSSBPLR(String tid1, String tid2, Map<String, Set<String>> Tid2goname,
			Map<String, Double> SSBPLR_map) throws Exception{
		int SSBP = 0;
		if (Tid2goname.get(tid1) != null && Tid2goname.get(tid2) != null){
			Set<String> gos_intersection = new HashSet<String>(Tid2goname.get(tid1));
			//Set<String> gos_i = Tid2goname.get(tid1);
			Set<String> gos_j = Tid2goname.get(tid2);
			gos_intersection.retainAll(gos_j);
			if (gos_intersection.size() != 0){
				int[] SSBP_array = new int[gos_intersection.size()+1];
				int go_count = 0;
				for (String go: gos_intersection) {
					String[] items = go.split(" ");
					go = items[0];
					if (m_goid2tid.get(go).size() > 0){
						go_count++;
						SSBP_array[go_count] = m_goid2tid.get(go).size();
					}
				}
				SSBP = getMin(SSBP_array);
			}
		}
		String SSBP_str = "0";
		if (SSBP == 0) SSBP_str = "1000000";
		else if (SSBP < 10) SSBP_str = "10";
		else if (SSBP < 50) SSBP_str = "50";
		else if (SSBP < 100) SSBP_str = "100";
		else if (SSBP < 500) SSBP_str = "500";
		else if (SSBP < 1000) SSBP_str = "1000";
		else if (SSBP < 1000000) SSBP_str = "1000000";
		double SSBP_LR = 1.0;
		if (SSBPLR_map.get(SSBP_str)!=null){
			SSBP_LR = SSBPLR_map.get(SSBP_str);
		}
		return SSBP_LR;
	}
	
	double getOrthologLR(String tid1, String tid2, Map<String, Double> RatioMap, 
			Map<String, Double> Ortholog_map) throws Exception {
		double Class = 0.0;
		if (HomologyTid_set.contains(tid1) && HomologyTid_set.contains(tid2)){
			String PPI = tid1 + ":" + tid2;
			if (RatioMap.get(PPI) != null)
				Class = RatioMap.get(PPI);
		}
		String Class_str = "0";
		if (Class > 0.8) Class_str = "0.8";
		else if (Class > 0.6) Class_str = "0.6";
		else if (Class > 0.4) Class_str = "0.4";
		else if (Class > 0.2) Class_str = "0.2";
		else if (Class > 0.0) Class_str = "0.0";
		double Ortholog_LR = 1.0;
		if (Ortholog_map.get(Class_str)!=null){
			Ortholog_LR = Ortholog_map.get(Class_str);
		}
		return Ortholog_LR;
	}
	
	void validation(List<ExpProfile> profile_list,
			List<Map<String, Double>> CoExpLR_map_list,
			Map<String, Double> DoEnrichLR_map, Map<String, Double> SSBPLR_map,
			Map<String, Double> OrthologLR_map) throws Exception {
		//int bin_size = 0;
		//String validation_fileName = null;
		/*for (String c: CoExpLR_map.keySet()){
			System.out.println(c +" "+CoExpLR_map.get(c));
		}*/
		for (String d: DoEnrichLR_map.keySet()){
			System.out.println(d +" "+DoEnrichLR_map.get(d));
		}
		for (String s: SSBPLR_map.keySet()){
			System.out.println(s +" "+SSBPLR_map.get(s));
		}
		for (String o: OrthologLR_map.keySet()){
			System.out.println(o +" "+OrthologLR_map.get(o));
		}
		//debug_println("Known2geneid " + m_Known2geneidsMap.size());
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(validation_option + "_interaction.log"));

		Map<String, Set<String>> tid2Pfam = new HashMap<String, Set<String>>();
		if (!PfamFile_fileName.equals("null")){
			tid2Pfam = getTid2PfamFile(PfamFile_fileName);
		}else{
			tid2Pfam = getTid2Pfam();
		}
		Map<String, Double> RatioMap = getRatioMap(RatioMap_fileName);
		Map<String, Set<String>> Tid2goname = new HashMap<String, Set<String>>();
		if (!PredictionGO_fileName.equals("null")){
			Tid2goname = getTid2gonamePrediction(PredictionGO_fileName);
		}else{
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0728_2010");
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0723_2009");
			Tid2goname = getTid2gonameUniProt("uniprot2go_0405_2011_BP");
		}	
		Map<String, Double> HomologyRatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);
		debug_println("sd_th "+ sd_th);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		for(int i=0; i < tid_list.size(); i++) {
			if (i % 1000 == 0)
				debug_println("#tid "+ i);
			String tid1 = tid_list.get(i);
			for(int j=0; j < tid_list.size(); j++) {
				String tid2 = tid_list.get(j);
				if (tid1.equals(tid2)){
					continue;
				}
				String tid2tid = tid1 + ":" + tid2;
				if (!Train_GSP_set.contains(tid2tid) 
						&& !Train_GSN_set.contains(tid2tid) 
						&& !Test_GSP_set.contains(tid2tid) 
						&& !Test_GSN_set.contains(tid2tid) 
								  ){
					continue;					
				}
				double LR = 0.0;
				if (validation_option.equals("all")){
					double CoExp_LR = 1.0;
					if (test_all_exp_flag){
						CoExp_LR = getExpMaxLR(i, j, profile_list, CoExpLR_map_list);
					}else{
						CoExp_LR = getExpLR(i, j, profile_list, CoExpLR_map_list);
					}
					double DoEnrich_LR = getDomainLR(tid1, tid2, tid2Pfam, RatioMap, DoEnrichLR_map);
					double SSBP_LR = getSSBPLR(tid1, tid2, Tid2goname, SSBPLR_map);
					double Ortholog_LR = getOrthologLR(tid1, tid2, HomologyRatioMap, OrthologLR_map);
					//System.out.println(" CoExp_LR" + CoExp_LR + " DoEnrich_LR" + DoEnrich_LR + " SSBP_LR" + SSBP_LR);
					LR = CoExp_LR * DoEnrich_LR * SSBP_LR * Ortholog_LR;
				}
				if (validation_option.equals("CoExp")){
					if (test_all_exp_flag){
						LR = getExpMaxLR(i, j, profile_list, CoExpLR_map_list);
					}else{
						LR = getExpLR(i, j, profile_list, CoExpLR_map_list);
					}
				}
				if (validation_option.equals("DoEnrich")){
					LR = getDomainLR(tid1, tid2, tid2Pfam, RatioMap, DoEnrichLR_map);
				}
				if (validation_option.equals("SSBP")){
					LR = getSSBPLR(tid1, tid2, Tid2goname, SSBPLR_map);
				}
				if (validation_option.equals("Ortholog")){
					LR = getOrthologLR(tid1, tid2, HomologyRatioMap, OrthologLR_map);
				}
				
				/*if (validation_option.equals("withoutOrtholog")){
					LR = CoExp_LR * DoEnrich_LR * SSBP_LR;
				}
				if (validation_option.equals("withoutSSBP_Ortholog")){
					LR = CoExp_LR * DoEnrich_LR;
				}
				if (validation_option.equals("random")){
					LR = random.nextInt(10000);
				}*/
	
				outInteractionFile(LR, tid1, tid2, bw1);
			}
		}
		bw1.close();
		//outValidation(gsp, gsn, test, total, bin_size, all_LR_map, validation_fileName);
	}

	static void outInteractionFile(double LR, String tid1, String tid2, BufferedWriter bw1) throws Exception{
		if (LR > LR_th){
			String geneid1 = m_Known2geneidsMap.get(tid1);
			bw1.write(tid1 + "\t");
			if (geneid1 != null) {
				//bw1.write(GeneAnnot.getGeneSym(geneid1) + "\t" + GeneAnnot.getGeneDesc(geneid1) + "\t");
				bw1.write(GeneAnnot.getGeneSym(geneid1) + "\t");
			}else{
				bw1.write("null\t");
			}
			String geneid2 = m_Known2geneidsMap.get(tid2);
			bw1.write(tid2 + "\t");
			if (geneid2 != null) {
				//bw1.write(GeneAnnot.getGeneSym(geneid2) + "\t" + GeneAnnot.getGeneDesc(geneid2) + "\t");
				bw1.write(GeneAnnot.getGeneSym(geneid2) + "\t");
			}else{
				bw1.write("null\t");
			}
			bw1.write(LR + "\n");
		}
	}

	static void getMaxMinLR(String fileName) throws Exception {
		Double LR = 0.0;
		max_LR = 0.0;
		min_LR = 1000000.0;	
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		while ((line = br.readLine()) != null) {
			String[] items = line.split("\t");
			if (items.length == 5) {
				LR = Double.parseDouble(items[4]);
				if (LR < min_LR) min_LR = LR;
				if (LR > max_LR) max_LR = LR;
			}
		}
		br.close();
	}	
	
	static void outPrecisionRecall(String fileName) throws Exception {
		//Set<String> single_TEST_set = getSingle_TESTset(TEST_set);
		String tid1 = null;
		String tid2 = null;
		Double LR = 0.0;
		if (auto_getMaxMinLR_flag) {
			getMaxMinLR(fileName);
			debug_println("auto set max_LR "+ max_LR + "  min_LR " + min_LR);
		}else{
			debug_println("manual set max_LR "+ max_LR + "  min_LR " + min_LR);
		}
		
		int bin_size = 1000;
		int[] prediction = new int[bin_size+1];
		int[] sel_prediction = new int[bin_size+1];
		int[] sel_GSN = new int[bin_size+1];
		int[] match = new int[bin_size+1];
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		while ((line = br.readLine()) != null) {
			String[] items = line.split("\t");
			if (items.length == 5) {
				tid1 = items[0];
				tid2 = items[2];
				LR = Double.parseDouble(items[4]);
				//int index = (int)(LR / 10);
				int index = (int)((LR - min_LR) / (max_LR - min_LR) * bin_size);
				if (index >= bin_size+1) 
					//continue;
					index = bin_size;
				boolean sel_prediction_found = false;
				boolean sel_GSN_found = false;
				boolean match_found = false;
				
				if ((Test_GSP_set.contains(tid1 + ":" + tid2))||(Test_GSP_set.contains(tid2 + ":" + tid1))){
					match_found = true;
				}
				
				if (single_Test_GSP_set.contains(tid1) && single_Test_GSP_set.contains(tid2)){
					sel_prediction_found = true;
					if ((Test_GSN_set.contains(tid1 + ":" + tid2))||(Test_GSN_set.contains(tid2 + ":" + tid1))){
						sel_GSN_found = true;
					}
				}
				
				for (int i = 0; i <= index; i++){
					prediction[i]++;
					if (sel_prediction_found) sel_prediction[i]++;
					if (match_found) match[i]++;
					if (sel_GSN_found) sel_GSN[i]++;
				}
			}
		}
		br.close();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("PrecisionRecall_bin"+bin_size+".log"));
		bw.write("index\tth\t#prediction\t#selected prediction\t#match\tprecision1\tprecision2\trecall\t-");
		bw.write("\n");
		
		for (int i = 0; i < prediction.length; i++){
			double th = (((double)i / bin_size) * (max_LR - min_LR) + min_LR );			
			double precision1 = (double) match[i] / (sel_GSN[i] + match[i]); 
			double precision2 = (double) match[i] / sel_prediction[i];
			double recall = (double) match[i] / Test_GSP_set.size();
			bw.write(i+"\t"+th+"\t"+prediction[i]+"\t"+sel_prediction[i]+"\t"+match[i]+"\t"+precision1+"\t"+precision2+"\t"+recall);
			bw.write("\n");
		}	
		bw.close();
		debug_println("Test_GSP_set.size() "+ Test_GSP_set.size());
	}
	
	Map<String, Set<String>> getPfam2tid() throws Exception {
		Map<String, Set<String>> Pfam2tid = new HashMap<String, Set<String>>();
		PreparedStatement getTid = m_conn
			.prepareStatement("SELECT distinct name, value"
				+ "	FROM ucsc_hg18.knownToPfam");		
		ResultSet rs = getTid.executeQuery();			
		while (rs.next()) {
			String tid = rs.getString("name");
			String Pfam = rs.getString("value");		
			Set<String> tid_set = Pfam2tid.get(Pfam);
			if (tid_set == null){
				tid_set = new HashSet<String>();
			}
			tid_set.add(tid);
			Pfam2tid.put(Pfam, tid_set);
		}		
		getTid.close();
		debug_println("getPfam2Tid() Pfam2Tid "+ Pfam2tid.size());
		return Pfam2tid;
	}
	
	static Map<String, Set<String>> getTid2Pfam() throws Exception {
		Map<String, Set<String>> tid2Pfam = new HashMap<String, Set<String>>();
		PreparedStatement getTid = m_conn
			.prepareStatement("SELECT distinct name, value"
				+ "	FROM ucsc_hg18.knownToPfam");		
		ResultSet rs = getTid.executeQuery();			
		while (rs.next()) {
			String tid = rs.getString("name");
			String Pfam = rs.getString("value");		
			Set<String> Pfam_set = tid2Pfam.get(tid);
			if (Pfam_set == null){
				Pfam_set = new HashSet<String>();
			}
			Pfam_set.add(Pfam);
			tid2Pfam.put(tid, Pfam_set);
		}		
		getTid.close();
		debug_println("getTid2Pfam() tid2Pfam "+ tid2Pfam.size());
		return tid2Pfam;
	}
	
	Map<String, Set<String>> getPfamFile2tid(String filename) throws Exception {
		debug_println("getPfamFile2tid() filename "+ filename);	
		Map<String, Set<String>> Pfam2tid = new HashMap<String, Set<String>>();
		List<String> Tid2Pfam = FileUtil.getFileAsList(filename);
		for (String line: Tid2Pfam) {
			String[] items = line.split("\t"); 
			String tid = items[0];		
			String Pfam = items[1];			
			Set<String> tid_set = Pfam2tid.get(Pfam);
			if (tid_set == null){
				tid_set = new HashSet<String>();
			}
			tid_set.add(tid);
			Pfam2tid.put(Pfam, tid_set);
		}		
		debug_println("getPfamFile2tid() PfamFile2Tid "+ Pfam2tid.size());
		return Pfam2tid;
	}
	
	static Map<String, Set<String>> getTid2PfamFile(String filename) throws Exception {
		debug_println("getTid2PfamFile() filename "+ filename);	
		Map<String, Set<String>> tid2Pfam = new HashMap<String, Set<String>>();		
		List<String> Tid2Pfam = FileUtil.getFileAsList(filename);
		for (String line: Tid2Pfam) {
			String[] items = line.split("\t"); 
			String tid = items[0];		
			String Pfam = items[1];		
			Set<String> Pfam_set = tid2Pfam.get(tid);
			if (Pfam_set == null){
				Pfam_set = new HashSet<String>();
			}
			Pfam_set.add(Pfam);
			tid2Pfam.put(tid, Pfam_set);
		}		
		debug_println("getTid2PfamFile() tid2PfamFile "+ tid2Pfam.size());
		return tid2Pfam;
	}
	
	static Map<String, Set<String>> getTid2gonamePrediction(String filename) 
	throws Exception {
		debug_println("getTid2gonamePrediction() filename "+ filename);		
		m_goid2tid = new HashMap<String, Set<String>>();
		Map<String, Set<String>> tid2goname = new HashMap<String, Set<String>>();
		List<String> prediction = FileUtil.getFileAsList(filename);
		for (String line: prediction) {
			String[] items = line.split("\t"); 
			String tid = items[1];
			//String uniprot = rs.getString("k.proteinID");			
			String goid = items[0];
			if (Enrich.getGoName(goid) == null) {
				continue;
			}
			{
				Set<String> go = tid2goname.get(tid);
				if (go == null) {
					go = new HashSet<String>();					
					tid2goname.put(tid, go);
				}
				String goname = Enrich.getGoName(goid);
				if (goname != null) {
					go.add(goname);
				}
			}
			{
				Set<String> tids = m_goid2tid.get(goid);
				if (tids == null) {
					tids = new HashSet<String>();
					m_goid2tid.put(goid, tids);
				}
				tids.add(tid);
			}			
		}
		debug_println("getTid2gonamePrediction() tid2goname "+ tid2goname.size());
		debug_println("  m_goid2tid "+ m_goid2tid.size());
		return tid2goname;
	}	
	
	static void debug_println(String str) {
		if (m_debug_flag) {
			System.out.println(str);
		}
	}

	static void debug_print(String str) {
		if (m_debug_flag) {
			System.out.print(str);
		}
	}
	
	Set<String> getRandomSet(Set<String> locs, int size) {
		List<String> list = new ArrayList<String>(locs);
		Set<String> result = new HashSet<String>();
		while (result.size() < size) {
			int i = random.nextInt(10000);
			result.add(list.get(i));
		}
		return result;
	}
}
