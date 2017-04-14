package edu.nchu.syslab.ppi;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.DriverManager;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.nchu.syslab.rnaseq.ExpProfile;
import edu.usc.zhoulab.common.File.FileUtil;
//import edu.usc.zhoulab.common.annot.GeneAnnot;

/**
 * 
 * @author ting
 * Prepare data:
 * 1. Training data (not one-to-one mapping): Label,Feature1,Feature2~8
 * 2. Test data (one-to-one mapping): Label,Feature1,Feature2~8
 * 
 * 1. output correlation(x,y) for both training and test sets
 *    correlation(x, y) means using set of label and set of one feature to calculate correlation.
 *    x is one of features (GO, RNA-seq, domain, or ortholog)
 *    y is the label: 0 no interaction, 1 interaction.
 * 2. perform logistic regression
 * 3. output correlation(f(x),y) for both training and test sets
 *    f(x) is the prediction from regression.
 * 
 * Please use the following configurations:
 * 1. one2one training and one2one test sets.
 * 2. all training and all test sets.
 * 3. all training and one2one test sets.
 * 
 * 
 */
public class OutputCorrData extends IsoIntActRawData{ 
	static double th_test = 0.1;
	static String Interaction_fileName = null;

	public static void main(String[] args) throws Exception {
		if (args.length != 8) {
			System.out.println("Usage: OutputCorrData option" +
					" Tid_fileName" +
					" RatioMap_fileName " +
					" PredictionGO_fileName " +
					" PfamFile_fileName" +
					" HomologyRatioMap_fileName" +
					" th_test GSPGSN_folderName");
			// execJava edu.nchu.syslab.ppi.SVMIsoIntAct PrecisionRecall\
			// /home/ting/IntAct/sel_tids.log \
			// /home/ting/IntAct/DomainRatioMap/RatioMap_file.log \
			// null \
			// null \
			// /home/jimliu/Annot/human2all/knownInt_t4_0.05.txt \
			// 0.0 true true > n_OutputCorrData.log
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
		//GSP_one2one_flag = Boolean.parseBoolean(args[7]);
		//TEST_one2one_flag = Boolean.parseBoolean(args[8]);
		GSPGSN_folderName = args[7];
		System.out.println("option " + option);
		System.out.println("Tid_fileName " + Tid_fileName);
		System.out.println("RatioMap_fileName " + RatioMap_fileName);
		System.out.println("PredictionGO_fileName " + PredictionGO_fileName);
		System.out.println("PfamFile_fileName " + PfamFile_fileName);
		System.out.println("HomologyRatioMap_fileName " + HomologyRatioMap_fileName);
		System.out.println("th_test " + th_test);
		//System.out.println("GSP_one2one_flag " + GSP_one2one_flag);
		//System.out.println("TEST_one2one_flag " + TEST_one2one_flag);
		System.out.println("GSPGSN_folderName " + GSPGSN_folderName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new OutputCorrData();
	}

	public OutputCorrData() throws Exception {
		//TODO: prepare data to run R
		prepare_all_arr(Profile_list);
		//RegressionIsoIntAct.ExcuteR();
		
		//TODO: get the coefficient result to do prediction and validation
		//List<Double> CoefList = RegressionIsoIntAct.getCoefList();
		
	}
	
	void prepare_all_arr(List<ExpProfile> profile_list) throws Exception {		
		//Set<String> single_old_GSP_set = getSingle_tid_set(Train_GSP_set);
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
		//Map<String, String> protein2knowns = GeneAnnot.getKnown2protein(hg18_to_uniprot);		
		debug_println("sd_th "+ sd_th);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		// create a test_list_isoforms
		List<String> datasets = FileUtil.getFileAsList("/home/jimliu/RNA-seq/Results/dataset.lst");
		String ds_str = datasets.toString().replace("[", "").replace("]", "").replace(" ", "");
		BufferedWriter bw_train = new BufferedWriter(new FileWriter("train_set.csv"));
		// d1, d2, d3, d4, d5, d7, d9, d10 and d11 
		//bw_train.write("Label,exp_d1_score,exp_d2_score,exp_d3_score,exp_d4_score,exp_d5_score,exp_d7_score,exp_d9_score,exp_d10_score,exp_d11_score,domain_score,GO_score,ortholog_score");
		bw_train.write("Label,"+ds_str+",domain_score,GO_score,ortholog_score");
		bw_train.write("\n");
		
		BufferedWriter bw_test = new BufferedWriter(new FileWriter("test_set.csv"));
		//bw_test.write("Label,exp_d1_score,exp_d2_score,exp_d3_score,exp_d4_score,exp_d5_score,exp_d7_score,exp_d9_score,exp_d10_score,exp_d11_score,domain_score,GO_score,ortholog_score");
		bw_test.write("Label,"+ds_str+",domain_score,GO_score,ortholog_score");
		bw_test.write("\n");
		
		int num_train_pos = 0;
		int num_train_neg = 0;
		//int num_train_both_pos_neg = 0;
		int num_test_pos = 0;
		int num_test_neg = 0;
		//int num_test_both_pos_neg = 0;
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
			
				String R_str = RegressionIsoIntAct.getExpScore_str(i, j, profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				//System.out.println(R + "	" + maxRatioD + "	" + SSBP + "	" + Class);
				
				if (Train_GSP_set.contains(tid2tid)){
					bw_train.write("1," + R_str +  maxRatioD + "," + SSBP + "," + Class); 
					bw_train.write("\n");
					num_train_pos++; // training positive samples
				}
				if (Train_GSN_set.contains(tid2tid)){
					//if (num_train_pos > num_train_neg) 
					//if (maxRatioD == 0)
					{
						bw_train.write("0," + R_str +  maxRatioD + "," + SSBP + "," + Class); 
						bw_train.write("\n");
						num_train_neg++; // training negative samples
					}
				}
				
				/*if (old_GSP_set.contains(tid2tid) && 
						((train_GSN_set.tids_membrane.contains(tid1) &&  train_GSN_set.tids_nucleus.contains(tid2)) ||
						(train_GSN_set.tids_membrane.contains(tid2) &&  train_GSN_set.tids_nucleus.contains(tid1)))){
					num_train_both_pos_neg++;
				}*/
				
				if (Test_GSP_set.contains(tid2tid)) {
					bw_test.write("1," + R_str + maxRatioD + "," + SSBP + "," + Class);
					bw_test.write("\n");
					num_test_pos++; // testing positive samples
				}
				if (Test_GSN_set.contains(tid2tid)) {
					//if (maxRatioD == 0) 
					{
						bw_test.write("0," + R_str + maxRatioD + "," + SSBP + "," + Class);
						bw_test.write("\n");
						num_test_neg++; // test negative samples
					}
				}
				
				/*if (TEST_set.contains(tid2tid) && (!TEST_set.contains(tid2tid)) &&
					(single_TEST_set.contains(tid1) && single_TEST_set.contains(tid2)
							&& ((GSN_set.tids_membrane.contains(tid1) && GSN_set.tids_nucleus.contains(tid2)) || 
								(GSN_set.tids_membrane.contains(tid2) && GSN_set.tids_nucleus.contains(tid1))) ) )
				{
					num_test_both_pos_neg++;
				}*/
			}
		}
		
		//debug_println("  train_one2one_flag " + GSP_one2one_flag);
		//debug_println("  test_one2one_flag " + TEST_one2one_flag);
		debug_println("  num_train_pos " + num_train_pos);
		debug_println("  num_train_neg " + num_train_neg);	
		//debug_println("  num_train_both_pos_neg " + num_train_both_pos_neg);
		debug_println("  num_test_pos " + num_test_pos);
		debug_println("  num_test_neg " + num_test_neg);
		//debug_println("  num_test_both_pos_neg " + num_test_both_pos_neg);
		bw_train.close();
		bw_test.close();
	}
	

	/*void out_all_arr(List<ExpProfile> profile_list, Set<String> knownInt_pos, Set<String> knownInt_neg) throws Exception {		
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
			Tid2goname = getTid2gonameUniProt("uniprot2go_0723_2009");
		}	
		Map<String, Double> HomologyRatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);
		Map<String, String> protein2knowns = GeneAnnot.getKnown2protein(hg18_to_uniprot);		
		debug_println("sd_th "+ sd_th);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
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
				if (pid1.contains("-"))
					pid1 = pid1.substring(0, pid1.indexOf('-'));
				if (pid2.contains("-"))
					pid2 = pid2.substring(0, pid2.indexOf('-'));
				String PPI = pid1 + ":" + pid2;
				if (!knownInt_pos.contains(tid2tid) && !knownInt_neg.contains(tid2tid)){
					continue;					
				}
			
				String R_str = RegressionIsoIntAct.getExpScore_str(i, j, profile_list);
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
	}*/
	
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

