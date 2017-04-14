package edu.nchu.syslab.ppi;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.usc.zhoulab.common.MathUtil;
import edu.nchu.syslab.rnaseq.ExpProfile;
import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.GeneAnnot;

/**
 * 
 * @author ting
 * Prepare data:
 * 1. Training data (not one-to-one mapping): Label \t ProteinID \t TranscriptID \t Features
 * 2. Test data (one-to-one mapping): Label \t ProteinID \t TranscriptID \t Features
 */
public class OutputMILdata extends IsoIntActRawData{ 
	static double th_test = 0.1;
	static String Interaction_fileName = null;

	public static void main(String[] args) throws Exception {
		if (args.length != 9) {
			System.out.println("Usage: OutputMILdata option" +
					" Tid_fileName" +
					" RatioMap_fileName " +
					" PredictionGO_fileName " +
					" PfamFile_fileName" +
					" HomologyRatioMap_fileName" +
					" th_test GSP_one2one_flag" +
					" Interaction_fileName" +
					" GSPGSN_folderName");
			// execJava edu.nchu.syslab.ppi.SVMIsoIntAct PrecisionRecall\
			// /home/ting/IntAct/sel_tids.log \
			// /home/ting/IntAct/DomainRatioMap/RatioMap_file.log \
			// null \
			// null \
			// /home/jimliu/Annot/human2all/knownInt_t4_0.05.txt \
			// 0.0 true > n_SVM.log
			// SetA
			return;
		}
		option = args[0];
		Tid_fileName = args[1];
		RatioMap_fileName = args[2];
		PredictionGO_fileName = args[3];
		PfamFile_fileName = args[4];
		HomologyRatioMap_fileName = args[5];
		th_test = Double.parseDouble(args[6]);
		GSP_one2one_flag = Boolean.parseBoolean(args[7]);
		GSPGSN_folderName = args[8];
		//Interaction_fileName = args[8];
		System.out.println("option " + option);
		System.out.println("Tid_fileName " + Tid_fileName);
		System.out.println("RatioMap_fileName " + RatioMap_fileName);
		System.out.println("PredictionGO_fileName " + PredictionGO_fileName);
		System.out.println("PfamFile_fileName " + PfamFile_fileName);
		System.out.println("HomologyRatioMap_fileName " + HomologyRatioMap_fileName);
		System.out.println("th_test " + th_test);
		System.out.println("GSP_one2one_flag " + GSP_one2one_flag);
		System.out.println("GSPGSN_folderName " + GSPGSN_folderName);
		//System.out.println("Interaction_fileName " + Interaction_fileName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new OutputMILdata();
		//new OutputMILdata("diverse");
	}

	public OutputMILdata() // Ting's program, the original version
			throws Exception {

		//List<ExpProfile> Profile_list = prepareExpProfile_list(tid_list);	
		if (option.equals("OutputMILdata")){
			prepare_all_arr();
		}
		if (option.equals("OutputMILdata2")){
			prepare_all_arr2();
		}
		if (option.equals("PrecisionRecall")){
			outPrecisionRecall(Interaction_fileName);
		}
		if (option.equals("Diverse")){
			//Set<String> all = getGSPsetIntActIsoJim(m_conn, "/home/ting/IntAct/2010_07_23/psimitab/intact.txt");
			IntSet Int_set = getGSPsetIntActIsoDiverse(m_conn, "/home/ting/IntAct/2011_03_09/psimitab/intact.txt");	
			Set<String> knownInt_pos = Int_set.knownInt_pos;
			Set<String> knownInt_neg = Int_set.knownInt_neg;
			out_all_arr(Profile_list, knownInt_pos, knownInt_neg);
		}
		if (option.equals("DiverseUseOldRule")){
			Set<String> knownInt_pos = getGSPsetIntActIsoJim(m_conn, "/home/ting/IntAct/2011_03_09/psimitab/intact.txt");
			out_all_arr_old(Profile_list, knownInt_pos);
		}
	}
	
	void out_all_arr(List<ExpProfile> profile_list, Set<String> knownInt_pos, Set<String> knownInt_neg) throws Exception {		
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
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
		Map<String, String> protein2knowns = GeneAnnot.getKnown2protein(hg18_to_uniprot);		
		debug_println("sd_th "+ sd_th);
		//List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		// create a test_list_isoforms
		
		BufferedWriter bw_train = new BufferedWriter(new FileWriter("IsoIntAct_pos_neg.out"));
		bw_train.write("Label\tProteinID\tTranscriptID\t" +
				"exp_d2_score\texp_d3_score\texp_d4_score\texp_d5_score\texp_d9_score\tdomain_score\tGO_score\tortholog_score");
		bw_train.write("\n");
		
		int num_train_pos = 0;
		int num_train_neg = 0;

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
				String pid1 = protein2knowns.get(tid1);
				String pid2 = protein2knowns.get(tid2);
				if (pid1 == null || pid2 == null){
					continue;
				}
				/*if (pid1.contains("-"))
					pid1 = pid1.substring(0, pid1.indexOf('-'));
				if (pid2.contains("-"))
					pid2 = pid2.substring(0, pid2.indexOf('-'));*/
				String PPI = pid1 + ":" + pid2;
				if (!knownInt_pos.contains(tid2tid) && !knownInt_neg.contains(tid2tid)){
					continue;					
				}
			
				String R_str = getExpScore_str(i, j, profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				//System.out.println(R + "	" + maxRatioD + "	" + SSBP + "	" + Class);
				
				if (knownInt_pos.contains(tid2tid)){	
					bw_train.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
					bw_train.write("\n");
					num_train_pos++; // training positive samples
				}
				if (knownInt_neg.contains(tid2tid)){	
					bw_train.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
					bw_train.write("\n");
					num_train_neg++; // training negative samples
				}
			}
		}
		
		debug_println("  num_train_pos " + num_train_pos);
		debug_println("  num_train_neg " + num_train_neg);	
		bw_train.close();
	}
	
	void out_all_arr_old(List<ExpProfile> profile_list, Set<String> knownInt_pos) throws Exception {		
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
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
			Tid2goname = getTid2gonameUniProt("uniprot2go_0405_2011");
		}	
		Map<String, Double> HomologyRatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);
		Map<String, String> protein2knowns = GeneAnnot.getKnown2protein(hg18_to_uniprot);		
		debug_println("sd_th "+ sd_th);
		Set<String> single_pos_set = getSingle_tid_set(knownInt_pos);
		
		BufferedWriter bw_train = new BufferedWriter(new FileWriter("IsoIntAct_pos_neg.out"));
		bw_train.write("Label\tProteinID\tTranscriptID\t" +
				"exp_d2_score\texp_d3_score\texp_d4_score\texp_d5_score\texp_d9_score\tdomain_score\tGO_score\tortholog_score");
		bw_train.write("\n");
		
		int num_train_pos = 0;
		int num_train_neg = 0;

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
				String pid1 = protein2knowns.get(tid1);
				String pid2 = protein2knowns.get(tid2);
				if (pid1 == null || pid2 == null){
					continue;
				}
				/*if (pid1.contains("-"))
					pid1 = pid1.substring(0, pid1.indexOf('-'));
				if (pid2.contains("-"))
					pid2 = pid2.substring(0, pid2.indexOf('-'));*/
				String PPI = pid1 + ":" + pid2;
				if (!knownInt_pos.contains(tid2tid) 
						&& !((GSN_set.tids_membrane.contains(tid1) && 
							  GSN_set.tids_nucleus.contains(tid2)) ||
							 (GSN_set.tids_membrane.contains(tid2) && 
							  GSN_set.tids_nucleus.contains(tid1)))
								  ){
					continue;					
				}
			
				String R_str = getExpScore_str(i, j, profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				//System.out.println(R + "	" + maxRatioD + "	" + SSBP + "	" + Class);
				
				if (knownInt_pos.contains(tid2tid)){	
					bw_train.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
					bw_train.write("\n");
					num_train_pos++; // training positive samples
				}else if (single_pos_set.contains(tid1) && single_pos_set.contains(tid2)){
					if ((GSN_set.tids_membrane.contains(tid1) && GSN_set.tids_nucleus.contains(tid2)) ||
						(GSN_set.tids_membrane.contains(tid2) && GSN_set.tids_nucleus.contains(tid1))){
						bw_train.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
						bw_train.write("\n");
						num_train_neg++; // training negative samples
					}
				}
			}
		}
		
		debug_println("  num_train_pos " + num_train_pos);
		debug_println("  num_train_neg " + num_train_neg);	
		bw_train.close();
	}
	
	void prepare_all_arr() throws Exception {		
		//Set<String> single_Train_GSP_set = getSingle_tid_set(Train_GSP_set);
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
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
		Map<String, String> protein2knowns = GeneAnnot.getKnown2protein(hg18_to_uniprot);		
		debug_println("sd_th "+ sd_th);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		// create a test_list_isoforms
		
		//GSN_set.tids_membrane.removeAll(single_Train_GSP_set);
		//GSN_set.tids_nucleus.removeAll(single_Train_GSP_set);
		
		BufferedWriter bw_train = new BufferedWriter(new FileWriter("IsoIntAct_train_set.out"));
		bw_train.write("Label\tProteinID\tTranscriptID\t" +
				"exp_d1_score\texp_d2_score\texp_d3_score\texp_d4_score\texp_d5_score\t" + 
				"exp_d7_score\texp_d9_score\texp_d10_score\texp_d11_score\t" + 
				"domain_score\tGO_score\tortholog_score");
		bw_train.write("\n");
		
		BufferedWriter bw_test = new BufferedWriter(new FileWriter("IsoIntAct_test_set.out"));
		bw_test.write("Label\tProteinID\tTranscriptID\t" +
				"exp_d1_score\texp_d2_score\texp_d3_score\texp_d4_score\texp_d5_score\t" + 
				"exp_d7_score\texp_d9_score\texp_d10_score\texp_d11_score\t" + 
				"domain_score\tGO_score\tortholog_score");
		bw_test.write("\n");
		
		int num_train_pos = 0;
		int num_train_neg = 0;
		int num_train_both_pos_neg = 0;
		int num_test_pos = 0;
		int num_test_neg = 0;
		int num_test_both_pos_neg = 0;
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
				String pid1 = protein2knowns.get(tid1);
				String pid2 = protein2knowns.get(tid2);
				if (pid1 == null || pid2 == null){
					continue;
				}
				/*if (pid1.contains("-"))
					pid1 = pid1.substring(0, pid1.indexOf('-'));
				if (pid2.contains("-"))
					pid2 = pid2.substring(0, pid2.indexOf('-'));*/
				String PPI = pid1 + ":" + pid2;
				if (!Train_GSP_set.contains(tid2tid) 
						&& !Train_GSN_set.contains(tid2tid) 
						&& !Test_GSP_set.contains(tid2tid) 
						&& !Test_GSN_set.contains(tid2tid) 
							  ){
					continue;					
				}
			
				String R_str = getExpScore_str(i, j, Profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				//System.out.println(R + "	" + maxRatioD + "	" + SSBP + "	" + Class);
				
				if (Train_GSP_set.contains(tid2tid)){	
					bw_train.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
					//bw_train.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + Class); 
					
					bw_train.write("\n");
					num_train_pos++; // training positive samples
				}
				if (Train_GSN_set.contains(tid2tid)){
					//if (num_train_pos > num_train_neg) 
					{
						bw_train.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
						//bw_train.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + Class); 
						bw_train.write("\n");
						num_train_neg++; // training negative samples
					}
				}
				
				if ((Train_GSP_set.contains(tid2tid)) && (Train_GSN_set.contains(tid2tid))){
					num_train_both_pos_neg++;
				}
				
				if (Test_GSP_set.contains(tid2tid)) {
					bw_test.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str + maxRatioD + "\t" + SSBP + "\t" + Class);
					//bw_test.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str + maxRatioD + "\t" + Class);
					bw_test.write("\n");
					num_test_pos++; // testing positive samples
				}
				if (Test_GSN_set.contains(tid2tid)) {
					bw_test.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str + maxRatioD + "\t" + SSBP + "\t" + Class);
					//bw_test.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str + maxRatioD + "\t" + Class);
					bw_test.write("\n");
					num_test_neg++; // test negative samples
				}
				
				if ((Test_GSP_set.contains(tid2tid)) && (Test_GSN_set.contains(tid2tid))){
					num_train_both_pos_neg++;
				}
			}
		}
		
		debug_println("  num_train_pos " + num_train_pos);
		debug_println("  num_train_neg " + num_train_neg);	
		debug_println("  num_train_both_pos_neg " + num_train_both_pos_neg);
		debug_println("  num_test_pos " + num_test_pos);
		debug_println("  num_test_neg " + num_test_neg);
		debug_println("  num_test_both_pos_neg " + num_test_both_pos_neg);
		bw_train.close();
		bw_test.close();
	}
	
	void prepare_all_arr2() throws Exception {		
		//Set<String> single_Train_GSP_set = getSingle_tid_set(Train_GSP_set);
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
		Map<String, Set<String>> pos_degree_map = new HashMap<String, Set<String>>();
		Map<String, Set<String>> neg_degree_map = new HashMap<String, Set<String>>();
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
		Map<String, String> protein2knowns = GeneAnnot.getKnown2protein(hg18_to_uniprot);		
		debug_println("sd_th "+ sd_th);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		// create a test_list_isoforms
		
		//GSN_set.tids_membrane.removeAll(single_Train_GSP_set);
		//GSN_set.tids_nucleus.removeAll(single_Train_GSP_set);
		
		BufferedWriter bw_pos = new BufferedWriter(new FileWriter("IsoIntAct_pos.out"));
		bw_pos.write("Label\tProteinID\tTranscriptID\t" +
				"exp_d1_score\texp_d2_score\texp_d3_score\texp_d4_score\texp_d5_score\t" + 
				"exp_d7_score\texp_d9_score\texp_d10_score\texp_d11_score\t" + 
				"domain_score\tGO_score\tortholog_score");
		bw_pos.write("\n");
		
		BufferedWriter bw_neg = new BufferedWriter(new FileWriter("IsoIntAct_neg.out"));
		bw_neg.write("Label\tProteinID\tTranscriptID\t" +
				"exp_d1_score\texp_d2_score\texp_d3_score\texp_d4_score\texp_d5_score\t" + 
				"exp_d7_score\texp_d9_score\texp_d10_score\texp_d11_score\t" + 
				"domain_score\tGO_score\tortholog_score");
		bw_neg.write("\n");
		
		int num_pos = 0;
		int num_neg = 0;
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
				String pid1 = protein2knowns.get(tid1);
				String pid2 = protein2knowns.get(tid2);
				if (pid1 == null || pid2 == null){
					continue;
				}
				/*if (pid1.contains("-"))
					pid1 = pid1.substring(0, pid1.indexOf('-'));
				if (pid2.contains("-"))
					pid2 = pid2.substring(0, pid2.indexOf('-'));*/
				String PPI = pid1 + ":" + pid2;
				if (!Train_GSP_set.contains(tid2tid) 
						&& !Train_GSN_set.contains(tid2tid) 
						&& !Test_GSP_set.contains(tid2tid) 
						&& !Test_GSN_set.contains(tid2tid) 
							  ){
					continue;					
				}

				if (Train_GSP_set.contains(tid2tid) || Test_GSP_set.contains(tid2tid)){	
					Set<String> tid2_set = pos_degree_map.get(tid1);
					if (tid2_set == null){
						tid2_set = new HashSet<String>();
					}
					tid2_set.add(tid2);
					pos_degree_map.put(tid1, tid2_set);
				}
				if (Train_GSN_set.contains(tid2tid) || Test_GSN_set.contains(tid2tid)){
					Set<String> tid2_set = neg_degree_map.get(tid1);
					if (tid2_set == null){
						tid2_set = new HashSet<String>();
					}
					tid2_set.add(tid2);
					neg_degree_map.put(tid1, tid2_set);
				}
			}
		}
		
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
				String pid1 = protein2knowns.get(tid1);
				String pid2 = protein2knowns.get(tid2);
				if (pid1 == null || pid2 == null){
					continue;
				}
				/*if (pid1.contains("-"))
					pid1 = pid1.substring(0, pid1.indexOf('-'));
				if (pid2.contains("-"))
					pid2 = pid2.substring(0, pid2.indexOf('-'));*/
				String PPI = pid1 + ":" + pid2;
				if (!Train_GSP_set.contains(tid2tid) 
						&& !Train_GSN_set.contains(tid2tid) 
						&& !Test_GSP_set.contains(tid2tid) 
						&& !Test_GSN_set.contains(tid2tid) 
							  ){
					continue;					
				}
			
				String R_str = getExpScore_str(i, j, Profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				//System.out.println(R + "	" + maxRatioD + "	" + SSBP + "	" + Class);
				
				if (Train_GSP_set.contains(tid2tid) || Test_GSP_set.contains(tid2tid)){	
					
					if (pos_degree_map.get(tid1).size()>=10){
						bw_pos.write("1\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 	
						bw_pos.write("\n");
						num_pos++;
					}

				}
				if (Train_GSN_set.contains(tid2tid) || Test_GSN_set.contains(tid2tid)){
			
					if (neg_degree_map.get(tid1).size()>=10){
						bw_neg.write("0\t" + PPI + "\t" + tid2tid + "\t" + R_str +  maxRatioD + "\t" + SSBP + "\t" + Class); 
						bw_neg.write("\n");
						num_neg++;
					}
				}
				
			}
		}
		
		debug_println("  num_pos " + num_pos);
		debug_println("  num_neg " + num_neg);	
		bw_pos.close();
		bw_neg.close();
	}
	
	static double getMaxExpScore(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		double R = 0;
		double maxCoExp_LR = -2;
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
			}
			else {
				R = 0;
			}
			if (maxCoExp_LR < R) {
				maxCoExp_LR = R;
			}
		}
		return maxCoExp_LR;
	}
	
	static String getExpScore_str_abs(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		double R = 0;
		String R_str = "";
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
			}
			else {
				R = 0;
			}
			R_str = R_str + R + "\t";
		}
		return R_str;
	}
	
	static String getExpScore_str(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		double R = 0;
		String R_str = "";
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = MathUtil.getPearson(exp[i], exp[j]);
			}
			else {
				R = 0;
			}
			R_str = R_str + R + "\t";
		}
		return R_str;
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

	static Set<String> getGSPsetIntActIsoJim(Connection conn, String fileName) throws Exception{
		// in protein level, the mapping is not one2one
		// in isoform level, the mapping is one2one
		Set<String> known_set = new HashSet<String>();
		Set<String> knownInt = new HashSet<String>();
		Set<String> knownInt_null = new HashSet<String>();
		Set<String> proteinInt = new HashSet<String>();
		//Map<String, String> known2protein = getKnown2protein(conn); 
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
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
			
			if (protein1.equals(protein2)) {
				continue;
			}

			if (!protein1.contains("-") || !protein2.contains("-")) {
				continue;
			}
			proteinInt.add(protein1+":"+protein2);	
			proteinInt.add(protein2+":"+protein1);
			
			all_proteins_set.add(protein1);
			all_proteins_set.add(protein2);
			
			Set<String> knowns1_one2one = protein2knowns.get(protein1);
			Set<String> knowns2_one2one = protein2knowns.get(protein2);
			
			debug_println("  protein1 "+ protein1 + "->"+knowns1_one2one);
			debug_println("  protein2 "+ protein2 + "->"+knowns2_one2one);
			String protein1org = protein1.substring(0, protein1.indexOf('-'));			
			String protein2org = protein2.substring(0, protein2.indexOf('-'));
			Set<String> knowns1 = protein2knowns.get(protein1org);
			Set<String> knowns2 = protein2knowns.get(protein2org);
			debug_println("    -> "+ protein1org +" -> " + knowns1);
			debug_println("    -> "+ protein2org +" -> " + knowns2);
			if (knowns1_one2one == null || knowns2_one2one == null) {
				knownInt_null.add(protein1 + ":" + protein2);
				knownInt_null.add(protein2 + ":" + protein1);				
				continue;
			}
			if (knowns1.size() > 1 && knowns2.size() > 1) {
				sel_proteins_set.add(protein1);
				sel_proteins_set.add(protein2);
				known_set.addAll(knowns1_one2one);
				known_set.addAll(knowns2_one2one);
				for (String g1 : knowns1_one2one) {
					for (String g2 : knowns2_one2one) {
						knownInt.add(g1 + ":" + g2);
						knownInt.add(g2 + ":" + g1);
					}
				}
			}
		}
		debug_println(fileName);
		//debug_println("one2one "+ one2one);
		debug_println("proteinInt "+ (proteinInt.size()/2)); //total human PPIs: 
		debug_println("all_proteins_set "+ all_proteins_set.size()); //human proteins: 
		debug_println("sel_proteins_set "+ sel_proteins_set.size()); //human proteins with one-to-one mapping to knownGene: 
		debug_println("known_set "+ (known_set.size())); //knownGenes in IntAct database:
		//FileUtil.outputCollectionAsFile(intActknowns, "sel_tids.log");
		debug_println("knownInt "+ (knownInt.size()/2)); //knownGene Interactions with one-to-one mapping to knownGene:
		debug_println("knownInt_null "+ (knownInt_null.size()/2)); 
		return knownInt;
	}
	
	/*public void OutputMILdataTing(String option) // Ting's program 
	throws Exception {
		debug_println("OutputMILdata() option " + option);
		Set<String> all = getGSPsetIntActIsoJim(m_conn, "/home/ting/IntAct/2010_07_23/psimitab/intact.txt");
		debug_println("all " + all.size());
	}*/
	
	class IntSet {
		Set<String> known_set = new HashSet<String>();
		Set<String> knownInt_pos = new HashSet<String>();
		Set<String> knownInt_neg = new HashSet<String>();
		Set<String> knownInt_null = new HashSet<String>();
		Set<String> proteinInt = new HashSet<String>();
		Set<String> all_proteins_set = new HashSet<String>();	
		Set<String> sel_proteins_set = new HashSet<String>();
	}
	
	IntSet getGSPsetIntActIsoDiverse(Connection conn, String fileName) throws Exception{
		// in protein level, the mapping is not one2one
		// in isoform level, the mapping is one2one

		//Map<String, String> known2protein = getKnown2protein(conn); 
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		List<String> lines = FileUtil.getFileAsList(fileName);
		IntSet Int_set = new IntSet();
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
			
			if (protein1.equals(protein2)) {
				continue;
			}

			if (!protein1.contains("-") || !protein2.contains("-")) {
				continue;
			}
			Int_set.proteinInt.add(protein1+":"+protein2);	
			Int_set.proteinInt.add(protein2+":"+protein1);
			
			Int_set.all_proteins_set.add(protein1);
			Int_set.all_proteins_set.add(protein2);
			
			Set<String> knowns1_one2one = protein2knowns.get(protein1);
			Set<String> knowns2_one2one = protein2knowns.get(protein2);
			
			debug_println("  protein1 "+ protein1 + "->"+knowns1_one2one);
			debug_println("  protein2 "+ protein2 + "->"+knowns2_one2one);
			String protein1org = protein1.substring(0, protein1.indexOf('-'));			
			String protein2org = protein2.substring(0, protein2.indexOf('-'));
			Set<String> knowns1 = protein2knowns.get(protein1org);
			Set<String> knowns2 = protein2knowns.get(protein2org);
			debug_println("    -> "+ protein1org +" -> " + knowns1);
			debug_println("    -> "+ protein2org +" -> " + knowns2);
			if (knowns1_one2one == null || knowns2_one2one == null) {
				Int_set.knownInt_null.add(protein1 + ":" + protein2);
				Int_set.knownInt_null.add(protein2 + ":" + protein1);				
				continue;
			}
			if (knowns1.size() > 1 && knowns2.size() > 1) {
				Int_set.sel_proteins_set.add(protein1);
				Int_set.sel_proteins_set.add(protein2);
				Int_set.known_set.addAll(knowns1_one2one);
				Int_set.known_set.addAll(knowns2_one2one);
				for (String g1 : knowns1_one2one) {
					for (String g2 : knowns2_one2one) {
						Int_set.knownInt_pos.add(g1 + ":" + g2);
						Int_set.knownInt_pos.add(g2 + ":" + g1);
					}
				}
				for (String g1 : knowns1){
					for (String g2 : knowns2){
						Int_set.knownInt_neg.add(g1 + ":" + g2);
						Int_set.knownInt_neg.add(g2 + ":" + g1);
					}
				}
			}
		}
		Int_set.knownInt_neg.removeAll(Int_set.knownInt_pos);
		debug_println(fileName);
		//debug_println("one2one "+ one2one);
		debug_println("proteinInt "+ (Int_set.proteinInt.size()/2)); //total human PPIs: 
		debug_println("all_proteins_set "+ Int_set.all_proteins_set.size()); //human proteins: 
		debug_println("sel_proteins_set "+ Int_set.sel_proteins_set.size()); //human proteins with one-to-one mapping to knownGene: 
		debug_println("known_set "+ (Int_set.known_set.size())); //knownGenes in IntAct database:
		//FileUtil.outputCollectionAsFile(intActknowns, "sel_tids.log");
		debug_println("knownInt_pos "+ (Int_set.knownInt_pos.size()/2)); //knownGene Interactions with one-to-one mapping to knownGene:
		debug_println("knownInt_neg "+ (Int_set.knownInt_neg.size()/2));
		debug_println("knownInt_null "+ (Int_set.knownInt_null.size()/2)); 
		return Int_set;
	}	
}
