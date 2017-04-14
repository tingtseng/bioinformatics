package edu.nchu.syslab.ppi;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.DriverManager;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.nchu.syslab.rnaseq.BuildExpMatrix;
import edu.nchu.syslab.rnaseq.ExpProfile;
import edu.usc.zhoulab.common.MathUtil;
import edu.usc.zhoulab.common.File.FileUtil;

public class IsoIntActRawData extends IsoIntAct{ 
	//static String m_pipeline = "wel_pipeline";
	//static String m_pipeline = "chao_pipeline";
	static String m_pipeline = "Jim_pipeline";
	static List<ExpProfile> Profile_list = null;
	
	public static void main(String[] args) 
	throws Exception
	{
		if (args.length != 12) {
			System.out.println("Usage: IsoIntAct option" +
					" folderName" +
					" Tid_fileName" +
					" RatioMap_fileName m_suffix" +
					" GSP_one2one_flag test_all_exp_flag Score_th" +
					" PredictionGO_fileName" +
					" PfamFile_fileName" +
					" HomologyRatioMap_fileName" +
					" validation_option" );
			//buildRatioMap
			// ENCODE_SRP000228_hg18_known_GM sel_tids.log buildRatioMap_n.log 0
			// true true 0.0
			// /home/jimliu/RNA-seq/isofrom_func/go2tid_p0.5_Oct2.txt 
			// /home/jimliu/RNA-seq/studyPfam/tid2pfam_f5.txt
			// /home/jimliu/RNA-seq/homology/knownInt_1e-6_0.1.txt
			// withoutSSBP_Ortholog
			return;
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
		System.out.println("option " + option);
		System.out.println("folderName " + folderName);
		System.out.println("Tid_fileName " + Tid_fileName);
		System.out.println("RatioMap_fileName " + RatioMap_fileName);
		System.out.println("m_suffix " + m_suffix);
		System.out.println("GSP_one2one_flag " + GSP_one2one_flag);
		System.out.println("test_all_exp_flag " + test_all_exp_flag);
		System.out.println("Score_th " + LR_th);
		System.out.println("PredictionGO_fileName " + PredictionGO_fileName);
		System.out.println("PfamFile_fileName " + PfamFile_fileName);
		System.out.println("HomologyRatioMap_fileName " + HomologyRatioMap_fileName);
		System.out.println("validation_option " + validation_option);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new IsoIntActRawData();
	}	
	
	public IsoIntActRawData() throws Exception {
		/*Set<String> old_GSP_set = getGSPsetIntActIso(m_conn, "/home/ting/IntAct/2009_07_31/psimitab/intact.txt", GSP_one2one_flag);
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true);
		m_Known2geneidsMap = GeneAnnot.getKnown2geneid(m_conn);
		m_Geneid2knownsMap = GeneAnnot.getGeneid2knowns(m_conn);
		GSNset GSN_set = getGSNsetIntActIso(m_Known2geneidsMap, m_Geneid2knownsMap);
		Set<String> TEST_set = getTESTset(); 
		//Set<String> single_TEST_set = getSingle_TESTset(TEST_set);
		//GSNset train_GSN_set = getGSNsetIntActIsoTraining(m_Known2geneidsMap, m_Geneid2knownsMap, single_TEST_set);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);*/
		Profile_list = prepareExpProfile_list(tid_list);
		//ExpProfile Profile = prepareExpProfile(tid_list);
		
		//TODO: option of validation (PrecisionRecall)
		if (option.equals("validationCoExp")) {
			if (test_all_exp_flag){
				System.out.println("test_all_exp_flag " + test_all_exp_flag);
				System.out.println("***use  Profile_list");
				validation(Train_GSP_set, GSN_set, Test_GSP_set, Profile_list);
				outPrecisionRecall(validation_option + "_interaction.log");
			}else{
				System.out.println("test_all_exp_flag " + test_all_exp_flag);
				System.out.println("***use  Profile");	
				//validation(old_GSP_set, GSN_set, TEST_set, Profile);
				outPrecisionRecall(validation_option + "_interaction.log");
			}
		}

		if (option.equals("validation")){
			validation(Train_GSP_set, GSN_set, Test_GSP_set, Profile_list);
			outPrecisionRecall(validation_option + "_interaction.log");
		}
		if (option.equals("PrecisionRecall")){
			outPrecisionRecall(validation_option + "_interaction.log");
		}
	}
	
	List<ExpProfile> prepareExpProfile_list(List<String> tid_list) throws Exception {		
		List<ExpProfile> Profile_list = new ArrayList<ExpProfile>();

		if (m_pipeline != null) {
			List<String> datasets = FileUtil
					.getFileAsList("/home/jimliu/RNA-seq/Results/" + m_pipeline
							+ "/dataset.lst");
			for (String d : datasets) {
				folderName = "/home/jimliu/RNA-seq/Results/" + m_pipeline + "/"
						+ d;
				debug_println("RNA-seq " + folderName);
				new BuildExpMatrix(folderName, m_conn);
				List<String> sampleList = getSampleList(folderName);
				Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix
						.getTransMap(sampleList);
				ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(
						sampleList, sample2expMap);
				debug_println("transExpProfile " + d + ": " + transExpProfile);
				transExpProfile.removeLowSD(sd_th); // sd_th 0.1
				transExpProfile.applyGeneList(tid_list);
				Profile_list.add(transExpProfile);
			}
		} else {
			for (int i = 1; i <= 11; i++) {
				if (i == 6 || i == 8)
					continue;
				folderName = "/home/jimliu/RNA-seq/Results/d" + i
						+ "_hg18_known";
				debug_println("RNA-seq " + folderName);
				new BuildExpMatrix(folderName, m_conn);
				List<String> sampleList = getSampleList(folderName);
				Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix
						.getTransMap(sampleList);
				ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(
						sampleList, sample2expMap);
				debug_println("transExpProfile d" + i + ": " + transExpProfile);
				transExpProfile.removeLowSD(sd_th); // sd_th 0.1
				transExpProfile.applyGeneList(tid_list);
				Profile_list.add(transExpProfile);
			}
		}
		debug_println("prepareExpProfile_list() " + Profile_list.size());
		return Profile_list;
	}
	
	ExpProfile prepareExpProfile(List<String> tid_list) throws Exception {
		int i = 3;
		folderName = "/home/jimliu/RNA-seq/Results/d"+i+"_hg18_known";
		new BuildExpMatrix(folderName, m_conn);
		List<String> sampleList = BuildExpMatrix.getSampleList(folderName);
		Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix.getTransMap(sampleList);
		ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(sampleList, sample2expMap);
		debug_println("transExpProfile d"+i+": " + transExpProfile);
		transExpProfile.removeLowSD(sd_th); // sd_th 0.1			
		transExpProfile.applyGeneList(tid_list);
		return transExpProfile;
	}
	
	/*ExpProfile prepareExpProfile(List<String> tid_list) throws Exception {
		new BuildExpMatrix(folderName, m_conn);
		List<String> sampleList = BuildExpMatrix.getSampleList(folderName);
		Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix.getTransMap(sampleList);
		ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(sampleList, sample2expMap);
		debug_println("transExpProfile: " + transExpProfile);
		transExpProfile.removeLowSD(sd_th); // sd_th 0.1			
		transExpProfile.applyGeneList(tid_list);
		return transExpProfile;
	}*/
	
	//TODO: validation of all data
	static double getMaxExpScore(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		double R = 0;
		double maxCoExp_Score = -2;
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = MathUtil.getPearson(exp[i], exp[j]);
			}
			else {
				R = 0;
			}
			if (maxCoExp_Score < R) {
				maxCoExp_Score = R;
			}
		}
		return maxCoExp_Score;
	}
	
	static double getExpScore(int i, int j, ExpProfile profile) 
	throws Exception {
		double R = 0;
			//int net = 1;//net = 0 -> co_exp_LR_d2; net = 1 -> co_exp_LR_d3
			double[][] exp = profile.exp;
			if (exp[i] != null && exp[j] != null) {
				R = MathUtil.getPearson(exp[i], exp[j]);
			}
		return R;
	}
	
	//finish
	static double getDomainScore(String tid1, String tid2, Map<String, Set<String>> tid2Pfam, 
			Map<String, Double> RatioMap) throws Exception {
		double maxRatioD = 0.0;
		double RatioD = 0.0;
		double do_int = 0.0;
		if (tid2Pfam.get(tid1) != null && tid2Pfam.get(tid2) != null){
			Set<String> domains_i = tid2Pfam.get(tid1);
			Set<String> domains_j = tid2Pfam.get(tid2);
			//double[] RatioD_array = new double[domains_i.size()*domains_j.size()+1];
			int do_count = 0;
			for (String di: domains_i) {
				for (String dj: domains_j){
					/*
					if (di.equals(dj))
						continue;*/
					do_count++;
					String domain_pair = di + ":" + dj;
					if (RatioMap.get(domain_pair) != null) {
						RatioD = RatioMap.get(domain_pair);
						do_int++;
					}
					if (maxRatioD < RatioD)
						maxRatioD = RatioD;
					if (maxRatioD > 100) maxRatioD = 100;
				}
			}
		}
		return maxRatioD;
		/*
		if (do_int > 0) {
			return 1;
		}
		return do_int; */		
	}

	//finish
	static double getGoScore(String tid1, String tid2, Map<String, Set<String>> Tid2goname) 
	throws Exception{
		double SSBP = 0;
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
				SSBP = 10.0 / (double)getMin(SSBP_array);
			}
		}
		return SSBP;
	}
	
	//finish
	static double getOrthologScore(String tid1, String tid2, Map<String, Double> RatioMap) throws Exception {
		double Class = 0.0;
		if (HomologyTid_set.contains(tid1) && HomologyTid_set.contains(tid2)){
			String PPI = tid1 + ":" + tid2;
			if (RatioMap.get(PPI) != null)
				Class = RatioMap.get(PPI);
		}
		return Class;
	}
	
	void validation(Set<String> old_GSP_set, GSNset GSN_set, Set<String> TEST_set, List<ExpProfile> profile_list)
	throws Exception { 
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
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
				String PPI = tid1 + ":" + tid2;
				if (!old_GSP_set.contains(PPI) 
						&& !TEST_set.contains(PPI) 
						&& !((GSN_set.tids_membrane.contains(tid1) && 
							  GSN_set.tids_nucleus.contains(tid2)) ||
							 (GSN_set.tids_membrane.contains(tid2) && 
							  GSN_set.tids_nucleus.contains(tid1)))
							  ){
					continue;					
				}
				
				double Score = 0.0;
				if (validation_option.equals("all")){
					double R = getMaxExpScore(i, j, profile_list);
					double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
					double SSBP = getGoScore(tid1, tid2, Tid2goname);
					double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
					//System.out.println(" CoExp_LR" + CoExp_LR + " DoEnrich_LR" + DoEnrich_LR + " SSBP_LR" + SSBP_LR);
					Score = R * maxRatioD * SSBP * Class;
				}
				if (validation_option.equals("CoExp")){
					Score = getMaxExpScore(i, j, profile_list);
				}
				if (validation_option.equals("DoEnrich")){
					Score = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				}
				if (validation_option.equals("SSBP")){
					Score = getGoScore(tid1, tid2, Tid2goname);
				}
				if (validation_option.equals("Ortholog")){
					Score = getOrthologScore(tid1, tid2, HomologyRatioMap);
				}
				
				/*if (validation_option.equals("withoutOrtholog")){
					Score = R * maxRatioD * SSBP;
				}
				if (validation_option.equals("withoutSSBP_Ortholog")){
					Score = R * maxRatioD;
				}
				if (validation_option.equals("random")){
					Score = random.nextInt(10000);
				}*/
				
				outInteractionFile(Score, tid1, tid2, bw1);
			}
		}
		bw1.close();
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
	
}
