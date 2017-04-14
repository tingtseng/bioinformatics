package edu.nchu.syslab.ppi;

import java.io.IOException;
import java.sql.DriverManager;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.usc.zhoulab.common.RocResult;
import edu.usc.zhoulab.common.File.FileUtil;

import edu.usc.zhoulab.common.MathUtil;
import edu.nchu.syslab.ppi.IsoIntActRawData;
import edu.nchu.syslab.rnaseq.ExpProfile;

import libsvm.svm_model;
import libsvm.svm_node;
import libsvm.svm_parameter;

public class SVMIsoIntAct extends IsoIntActRawData{ 
	static double th_test = 0.1;
	static String SVM_type = null;
	static int SVM_TYPE = 0;

	public static void main(String[] args) throws Exception {
		if (args.length != 9) {
			System.out.println("Usage: SVMIsoIntAct option" +
					" Tid_fileName" +
					" RatioMap_fileName " +
					" PredictionGO_fileName " +
					" PfamFile_fileName" +
					" HomologyRatioMap_fileName" +
					" th_test GSP_one2one_flag" +
					" SVM_TYPE");
			// execJava edu.nchu.syslab.ppi.SVMIsoIntAct null \
			// /home/ting/IntAct/sel_tids.log \
			// /home/ting/IntAct/DomainRatioMap/RatioMap_file.log \
			// null \
			// null \
			// /home/jimliu/Annot/human2all/knownInt_t4_0.05.txt \
			// 0.0 true > n_SVM.log
			// svm_parameter.RBF
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
		SVM_type = args[8];
		if (SVM_type.equals("POLY"))
			SVM_TYPE = svm_parameter.POLY;
		if (SVM_type.equals("LINEAR"))
			SVM_TYPE = svm_parameter.LINEAR;
		if (SVM_type.equals("RBF"))
			SVM_TYPE = svm_parameter.RBF;
		System.out.println("option " + option);
		System.out.println("Tid_fileName " + Tid_fileName);
		System.out.println("RatioMap_fileName " + RatioMap_fileName);
		System.out.println("PredictionGO_fileName " + PredictionGO_fileName);
		System.out.println("PfamFile_fileName " + PfamFile_fileName);
		System.out.println("HomologyRatioMap_fileName " + HomologyRatioMap_fileName);
		System.out.println("th_test " + th_test);
		System.out.println("GSP_one2one_flag " + GSP_one2one_flag);
		System.out.println("SVM_type " + SVM_type);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new SVMIsoIntAct();
	}

	public SVMIsoIntAct() throws Exception {
		/*Set<String> old_GSP_set = getGSPsetIntActIso(m_conn, "/home/ting/IntAct/2009_07_31/psimitab/intact.txt", GSP_one2one_flag);
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true);
		m_Known2geneidsMap = GeneAnnot.getKnown2geneid(m_conn);
		m_Geneid2knownsMap = GeneAnnot.getGeneid2knowns(m_conn);
		GSNset GSN_set = getGSNsetIntActIso(m_Known2geneidsMap, m_Geneid2knownsMap);
		Set<String> TEST_set = getTESTset(); 
		Set<String> single_TEST_set = getSingle_TESTset(TEST_set);
		GSNset train_GSN_set = getGSNsetIntActIsoTraining(m_Known2geneidsMap, m_Geneid2knownsMap, single_TEST_set);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		List<ExpProfile> Profile_list = prepareExpProfile_list(tid_list);*/
		
		Arr all_arr = prepare_all_arr();
		buildModel_SVM_byArrays2(all_arr.ptest_arr, all_arr.ntest_arr, all_arr.pos_arr, all_arr.neg_arr);
		
		//RocResult results = SVMPredictor2.svmPredict(model, nodes, expected);
		
		//results.writeToFile("/home/ting/IntAct/SVM/SVM_result.txt");
		
	}
	
	class Arr {
		ArrayList<double[]> pos_arr = new ArrayList<double[]>();
		ArrayList<double[]> neg_arr = new ArrayList<double[]>();
		ArrayList<double[]> ptest_arr = new ArrayList<double[]>();
		ArrayList<double[]> ntest_arr = new ArrayList<double[]>();
	}
	
	Arr prepare_all_arr()throws Exception{
		Arr all_arr = new Arr();
		
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
		debug_println("sd_th "+ sd_th);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		// create a test_list_isoforms
		
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
				if (!Train_GSP_set.contains(tid2tid) 
						&& !((GSN_set.tids_membrane.contains(tid1) && 
							  GSN_set.tids_nucleus.contains(tid2)) ||
							 (GSN_set.tids_membrane.contains(tid2) && 
							  GSN_set.tids_nucleus.contains(tid1)))
						&& !Test_GSP_set.contains(tid2tid) 
						&& !((GSN_set.tids_membrane.contains(tid1) && 
							  GSN_set.tids_nucleus.contains(tid2)) ||
							 (GSN_set.tids_membrane.contains(tid2) && 
							  GSN_set.tids_nucleus.contains(tid1)))
							  ){
					continue;					
				}
			
				List<Double> ExpScore_list = getExpScore_list(i, j, Profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				//System.out.println(R + "	" + maxRatioD + "	" + SSBP + "	" + Class);
				double[] e = new double[8];
				e[0] = ExpScore_list.get(0);
				e[1] = ExpScore_list.get(1);
				e[2] = ExpScore_list.get(2);
				e[3] = ExpScore_list.get(3);
				e[4] = ExpScore_list.get(4);
				e[5] = maxRatioD;
				e[6] = SSBP;
				e[7] = Class;
				
				if (Train_GSP_set.contains(tid2tid)){		
					all_arr.pos_arr.add(e); // training positive samples
				}
				if ((GSN_set.tids_membrane.contains(tid1) &&  GSN_set.tids_nucleus.contains(tid2)) ||
					(GSN_set.tids_membrane.contains(tid2) &&  GSN_set.tids_nucleus.contains(tid1))){
					if (all_arr.pos_arr.size() > all_arr.neg_arr.size()) {
						all_arr.neg_arr.add(e); // training negative samples
					}
				}
				
				if (Test_GSP_set.contains(tid2tid)){
					all_arr.ptest_arr.add(e); // testing positive samples
				}else if (single_Test_GSP_set.contains(tid1) && single_Test_GSP_set.contains(tid2) &&
						((GSN_set.tids_membrane.contains(tid1) &&  GSN_set.tids_nucleus.contains(tid2)) ||
						 (GSN_set.tids_membrane.contains(tid2) &&  GSN_set.tids_nucleus.contains(tid1))) ) {
					   all_arr.ntest_arr.add(e); // test negative samples
				}
			}
		}
		
		debug_println("  num_pos " + num_pos);
		debug_println("  num_neg " + num_neg);	
		debug_println("  all_arr.pos_arr_rank " + all_arr.pos_arr.size());
		debug_println("  all_arr.neg_arr_rank " + all_arr.neg_arr.size());		
		debug_println("  all_arr.ntest_arr_rank " + all_arr.ntest_arr.size());
		debug_println("  all_arr.ptest_arr_rank " + all_arr.ptest_arr.size());
		return all_arr;
	}
	
	static List<Double> getExpScore_list(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		double R = 0;
		List<Double> ExpScore_list = new ArrayList<Double>();
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
			}
			ExpScore_list.add(R);
		}
		return ExpScore_list;
	}
	
	void buildModel_SVM_byArrays2(ArrayList<double[]> ptest_arr_rank, ArrayList<double[]> ntest_arr_rank,
			ArrayList<double[]> pos_arr_rank, ArrayList<double[]> neg_arr_rank)
			throws IOException {

		List<svm_node[]> pnodes = get_svm_nodes2(pos_arr_rank);
		debug_println("  pnodes " + pnodes.size());
		List<svm_node[]> nnodes = get_svm_nodes2(neg_arr_rank);
		debug_println("  nnodes " + nnodes.size());

		svm_model model = SVMModelBuilder2.buildModel(pnodes, nnodes, SVM_TYPE, Math.pow(2,0), Math.pow(2,0));
		debug_println("  buildModel " + SVM_type + " " + SVM_TYPE + " C=" + Math.pow(2,0) + " gamma=" + Math.pow(2,0));
		/*{
			RocResult total = new RocResult(); // i.e. edu.usc.zhoulab.common.RocResult
			List<svm_node[]> test_nodes;
			test_nodes = get_svm_nodes2(ptest_arr_rank);
			RocResult x = SVMPredictor2.svmPredict(model, test_nodes, true);
			total.add(x);
			test_nodes = get_svm_nodes2(ntest_arr_rank);
			x = SVMPredictor2.svmPredict(model, test_nodes, false);
			total.add(x);
		}*/
		{
			List<svm_node[]> test_nodes;
			test_nodes = get_svm_nodes2(ptest_arr_rank);
			SVMPredictor2.svmPredict(model, test_nodes, true, 0);
			List<Double> ptest_arr_mean = new ArrayList<Double>();
			for (double[] arr: ptest_arr_rank) {
				ptest_arr_mean.add(MathUtil.getAver(arr));
			}
			Collections.sort(ptest_arr_mean);
			debug_println("ptest_arr_mean: "+ ptest_arr_mean);
		}
		
		/*double min_dec = SVMPredictor2.m_min_dec; 
		double max_dec = SVMPredictor2.m_max_dec;*/
		
		Set<Double> dec_set = new HashSet<Double>(SVMPredictor2.m_dec_list);
		List<Double> dec_list = new ArrayList<Double>(dec_set);
		Collections.sort(dec_list);
		//for(double th = min_dec; th <= max_dec; th += (max_dec - min_dec)/100)
		
		List<svm_node[]> test_nodes_pos = get_svm_nodes2(ptest_arr_rank, th_test);
		List<svm_node[]> test_nodes_neg = get_svm_nodes2(ntest_arr_rank, th_test);		
		debug_println("th_test "+ th_test +" test_nodes_pos " +test_nodes_pos.size());
		debug_println("th_test "+ th_test +" test_nodes_neg " +test_nodes_neg.size());
		
		debug_println("REP>\tth\tTP/TN/FP/FN\tTrueAcc\tRecall\tPrecision");
		for(double th: dec_list)
		{
			//double th = th1/1000;
			RocResult total = new RocResult(); // i.e. edu.usc.zhoulab.common.RocResult
			RocResult x = SVMPredictor2.svmPredict(model, test_nodes_pos, true, th);
			total.add(x);
			x = SVMPredictor2.svmPredict(model, test_nodes_neg, false, th);
			total.add(x);
			debug_println("REP>\t"+ th +"\t" +total.toStringDetail());
			// "TP/TN/FP/FN \t TrueAcc \t Recall \t Precision "
		}
		//return total;
	}

	List<svm_node[]> get_svm_nodes2(ArrayList<double[]> arrays) {
		List<svm_node[]> svm_nodes = new ArrayList<svm_node[]>();
		for (double[] arr : arrays) {
			svm_node[] nodes = new svm_node[arr.length];
			for (int j = 0; j < arr.length; j++) {
				svm_node node = new svm_node();
				node.index = j + 1;
				node.value = (double) arr[j];
				nodes[j] = node;
			}
			svm_nodes.add(nodes);
		}
		return svm_nodes;
	}

	List<svm_node[]> get_svm_nodes2(ArrayList<double[]> arrays, double th) {
		List<svm_node[]> svm_nodes = new ArrayList<svm_node[]>();
		for (double[] arr : arrays) {
			svm_node[] nodes = new svm_node[arr.length];
			for (int j = 0; j < arr.length; j++) {
				svm_node node = new svm_node();
				node.index = j + 1;
				node.value = (double) arr[j];
				nodes[j] = node;
			}
			if (MathUtil.getAver(arr) >= th) {
				svm_nodes.add(nodes);
			}
		}
		return svm_nodes;
	}
	
	/*public static RocResult svmPredict(svm_model model,
			List<svm_node[]> nodes, boolean expected) {
		RocResult results = new RocResult();
		for (svm_node[] n : nodes) {
			// double svm_predict = svm.svm_predict(model, n);
			double[] dec_values = new double[1];
			svm.svm_predict_values(model, n, dec_values);
			int label = (dec_values[0] > 0.0) ? 1 : -1;
			results.addSample(expected, label, dec_values[0]);
		}
		return results;
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
