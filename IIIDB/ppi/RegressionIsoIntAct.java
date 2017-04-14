package edu.nchu.syslab.ppi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.sql.DriverManager;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.GeneAnnot;

import edu.usc.zhoulab.common.MathUtil;
import edu.nchu.syslab.ppi.IsoIntActRawData;
import edu.nchu.syslab.rnaseq.ExpProfile;



public class RegressionIsoIntAct extends IsoIntActRawData{
	
	static Map<String, String> m_tidInt2pidInt = null;
	//static Map<String, String> m_tidInt2pidInt2 = null;
	static Set<String> m_all_knowns_set = null;
	static Set<String> m_pidInt_set = null;
	//static Set<String> m_pidInt_set2 = null;
	static Set<String> m_all_proteins_set = null;
	static Set<String> m_sel_proteins_set = null;
	
	public static void main(String[] args) throws Exception {
		if (args.length != 13) {
			System.out.println("Usage: IsoIntAct option" +
					" folderName" +
					" Tid_fileName" +
					" RatioMap_fileName m_suffix" +
					" GSP_one2one_flag test_all_exp_flag LR_th" +
					" PredictionGO_fileName" +
					" PfamFile_fileName" +
					" HomologyRatioMap_fileName" +
					" validation_option" +
					" GSPGSN_folderName");
			// sel_tids.log buildRatioMap_n.log
			// /home/jimliu/RNA-seq/isofrom_func/go2tid_p0.5_Oct2.txt 
			// /home/jimliu/RNA-seq/studyPfam/tid2pfam_f5.txt
			// /home/jimliu/RNA-seq/homology/knownInt_1e-6_0.1.txt
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
		GSPGSN_folderName = args[12];
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
		System.out.println("GSPGSN_folderName " + GSPGSN_folderName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new RegressionIsoIntAct();
	}

	public RegressionIsoIntAct() throws Exception {
		
		if (option.equals("Regression")) {
			// TODO: prepare data to run R
			prepare_all_arr();
			ExcuteR();

			// TODO: get the coefficient result to do prediction and validation
			List<Double> CoefList = getCoefList();
			validation(CoefList);
			//out_validation(CoefList);
			outPrecisionRecall("interaction.log");
		}
		
		if (option.equals("outPrecisionRecall")){
			List<Double> CoefList = getCoefList();
			validation(CoefList);
			outPrecisionRecall("interaction.log");
		} 
		
		if (option.equals("OutValidation")){
			prepare_all_arr();
			ExcuteR();
			List<Double> CoefList = getCoefList();
			out_validation(CoefList);
		}
		
		if (option.equals("ProcessComplexData")){
			ProcessComplexData();
		}
		
		if (option.equals("OutputComplexData")){
			outputComplexData();
		}
		
		if (option.equals("ComplexAnalysisRecallIII")){
			//prepare_all_arr();
			//ExcuteR();
			List<Double> CoefList = getCoefList();
			ComplexAnalysis(CoefList);
		}
		if (option.equals("ComplexAnalysisRandom")){
			ComplexAnalysis_random();
		}
		
		if (option.equals("ComplexAnalysisIII")){
			ComplexAnalysis_isoform();
		}
		if (option.equals("ComplexAnalysisPPI")){
			ComplexAnalysis_protein();
		}
		
		if (option.equals("ComplexAnalysisRandomIso")){
			ComplexAnalysis_random_isoform();
		}
		if (option.equals("ComplexAnalysisRandomPro")){
			ComplexAnalysis_random_protein();
		}
		if (option.equals("Output_sif_file")){
			outputNetwork3("/home/ting/IntAct/Regression/test_final_2.0/out_interaction_features.log");
			//outputNetwork();
			//outputNetwork2();
		}
	}

	void outputNetwork() throws Exception {
		//uc004fqi.1      6192    RPS4Y1  uc004anf.2      3190    HNRNPK  3.585944321288961       null
		//BufferedWriter bw_iii = new BufferedWriter(new FileWriter("III.sif"));
		//BufferedWriter bw_ppi = new BufferedWriter(new FileWriter("PPI.sif"));
		Set<String> iii_set = new HashSet<String>();
		Set<String> ppi_set = new HashSet<String>();
 		List<String> lines = FileUtil.getFileAsList("out_interaction_features.log");
		for (String s: lines){
			String[] items = s.split("\t"); 
			//bw_iii.write(items[0] + " pp " + items[3] + "\n");
			String tid1 = items[0];
			String tid2 = items[3];
			String pid1 = items[2];
			String pid2 = items[5];
			if (iii_set.contains(tid2 + " pp " + tid1))
				continue;
			if (ppi_set.contains(pid2 + " pp " + pid1))
				continue;
			String geneid1 = m_Known2geneidsMap.get(tid1);
			String geneid2 = m_Known2geneidsMap.get(tid2);
			if (geneid1 != null && geneid2 != null){
				if (geneid1.equals(geneid2))
					continue;
			}
			iii_set.add(tid1 + " pp " + tid2);
			ppi_set.add(pid1 + " pp " + pid2);
		}
		//bw_iii.close();
		FileUtil.outputCollectionAsFile(iii_set, "III.sif");
		FileUtil.outputCollectionAsFile(ppi_set, "PPI.sif");
	}
	
	void outputNetwork2() throws Exception {
		Set<String> GSP_IntAct_set = getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		Set<String> iii_set = new HashSet<String>();
		for (String s :GSP_IntAct_set){
			String[] items = s.split(":"); 
			String tid1 = items[0];
			String tid2 = items[1];
			if (iii_set.contains(tid2 + " pp " + tid1))
				continue;
			String geneid1 = m_Known2geneidsMap.get(tid1);
			String geneid2 = m_Known2geneidsMap.get(tid2);
			if (geneid1 != null && geneid2 != null){
				if (geneid1.equals(geneid2))
					continue;
			}
			iii_set.add(tid1 + " pp " + tid2);
		}
		FileUtil.outputCollectionAsFile(iii_set, "III.sif");
		
		Set<String> ppi_set = new HashSet<String>();
		for (String p :m_pidInt_set){
			String[] items = p.split(":");
			String pid1 = items[0];
			String pid2 = items[1];
			if (ppi_set.contains(pid2 + " pp " + pid1))
				continue;
			ppi_set.add(pid1 + " pp " + pid2);
		}
		FileUtil.outputCollectionAsFile(ppi_set, "PPI.sif");
	}
	
	// Format: GeneSym1(KnownGeneID1) pp GeneSym2(KnownGeneID2)
	void outputNetwork3(String fileName) throws Exception {
		debug_println("outputNetwork3() " + fileName);
		List<String> lines = FileUtil.getFileAsList(fileName);
		BufferedWriter bw = new BufferedWriter(new FileWriter("interaction_25.sif"));
		int i = 0;
		for (String line: lines) {
			String[] items = line.split("\t");
			String tid1 = items[0];
			String genesym1 = items[2];
			String tid2 = items[4].trim();
			String genesym2 = items[6];
			double score = Double.parseDouble(items[8]);
			if (score >= LR_th){
				bw.write(genesym1 + "(" + tid1 + ") pp " + genesym2 + "(" + tid2 + ")\n");
				i++;
				if (i<=50){
					debug_print(genesym1 + "(" + tid1 + ") pp " + genesym2 + "(" + tid2 + ")\n");
				}
			}
		}
		bw.close();
	}
	
	/**
	 * Protein complex analysis 
	 * 1. Process protein complex database 
	 * 2. Build 3 networks 
	 *    III networks with threshold 2.3710 (At recall 5%, the precision is 64.1%) 
     *    III networks with threshold 1.6970 (At recall 10%, the precision is 48.8%) 
     *    III networks with threshold 1.3479 (At recall 20%, the precision is 43.8%) 
     *    III networks with threshold 1.0394 (At recall 30%, the precision is 35.2%)
	 *    IntAct PPI.
	 * 3. Calculate (# interactions within complex) / (# all interactions)
	 */
	
	static Set<String> getGSPsetIntActIso(String fileName, boolean one2one) throws Exception{
		debug_println("getGSPsetIntActIso()");
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		m_tidInt2pidInt = new HashMap<String, String>();
		m_all_knowns_set = new HashSet<String>();
		m_pidInt_set = new HashSet<String>();
		m_all_proteins_set = new HashSet<String>();	
		m_sel_proteins_set = new HashSet<String>();
		List<String> lines = FileUtil.getFileAsList(fileName);

		for (String line: lines) {
			if (line.startsWith("#") || !line.contains("taxid:9606(Human)")) {
				continue;
			}
			String[] items = line.split("\t");
			String protein1 = items[0];
			String protein2 = items[1];
			if (!protein1.startsWith("uniprotkb:")) {
				continue;
			}
			if (!protein2.startsWith("uniprotkb:")) {
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
			
			m_all_proteins_set.add(protein1);
			m_all_proteins_set.add(protein2);
			
			Set<String> knowns1 = protein2knowns.get(protein1);
			Set<String> knowns2 = protein2knowns.get(protein2);
			//debug_println("  protein1 "+ protein1 +" -> knowns "+ knowns1);
			if (knowns1 == null && protein1.contains("-")) {
			//if (protein1.contains("-")) {
				protein1 = protein1.substring(0, protein1.indexOf('-'));
				knowns1 = protein2knowns.get(protein1);
				//debug_println("    "+ protein1 +" -> knowns "+ knowns1);
			}
			if (knowns1 == null) {
				continue;
			}
			if (knowns2 == null && protein2.contains("-")) {
			//if (protein2.contains("-")) {		
				protein2 = protein2.substring(0, protein2.indexOf('-'));
				knowns2 = protein2knowns.get(protein2);
			}
			if (knowns2 == null) {
				continue;
			}
			
			if (one2one) {
				if (knowns1.size() == 1 && knowns2.size() == 1) { // check one-to-one
					m_sel_proteins_set.add(protein1);
					m_sel_proteins_set.add(protein2);
					m_all_knowns_set.addAll(knowns1);
					m_all_knowns_set.addAll(knowns2);
					for (String g1 : knowns1) {
						for (String g2 : knowns2) {
							m_tidInt2pidInt.put(g1 + ":" + g2, protein1 +":" + protein2);
							m_tidInt2pidInt.put(g2 + ":" + g1, protein2 +":" + protein1);
							m_pidInt_set.add(protein1 +":" + protein2);
							m_pidInt_set.add(protein2 +":" + protein1);
						}
					}
				}
			} else {
				if (knowns1.size() >= 1 && knowns2.size() >= 1) {
					m_sel_proteins_set.add(protein1);
					m_sel_proteins_set.add(protein2);
					m_all_knowns_set.addAll(knowns1);
					m_all_knowns_set.addAll(knowns2);
					for (String g1 : knowns1) {
						for (String g2 : knowns2) {
							// build the map: tid-pair -> protein-pair
							m_tidInt2pidInt.put(g1 + ":" + g2, protein1 +":" + protein2);
							m_tidInt2pidInt.put(g2 + ":" + g1, protein2 +":" + protein1);
							m_pidInt_set.add(protein1 +":" + protein2);
							m_pidInt_set.add(protein2 +":" + protein1);
						}
					}
				}
			}
		}
		debug_println(fileName);
		debug_println("one2one "+ one2one);
		debug_println("proteinInt.size()/2 "+ (m_pidInt_set.size()/2)); //total human PPIs: 
		debug_println("all_proteins_set "+ m_all_proteins_set.size()); //human proteins: 
		debug_println("sel_proteins_set "+ m_sel_proteins_set.size()); //human proteins with one-to-one mapping to knownGene: 
		debug_println("knownInt.size()/2 "+ (m_tidInt2pidInt.keySet().size()/2)); //knownGene Interactions with one-to-one mapping to knownGene:
		debug_println("known_set "+ (m_all_knowns_set.size())); //knownGenes in IntAct database:
		return m_tidInt2pidInt.keySet();
	}
	
	void ComplexAnalysis_isoform() throws Exception {
		Set<String> complex_tid2tid_set = FileUtil.getFileAsSet("../complex_tid2tid_set_" + GSP_one2one_flag);
		Set<String> complex_tid_set = new HashSet<String>();
		for (String s: complex_tid2tid_set){
			String[] items = s.split(":"); 
			complex_tid_set.add(items[0]);
			complex_tid_set.add(items[1]);
		}
		debug_println(" # complex cover isoforms  = " + complex_tid_set.size() );
		Set<String> GSP_IntAct_set = getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
		int num_within_complex = 0;
		int num_withComplex = 0;
		for (String tid2tid: GSP_IntAct_set){
			String[] items = tid2tid.split(":");
			String tid1 = items[0];
			String tid2 = items[1];
			String geneid1 = m_Known2geneidsMap.get(tid1);
			String geneid2 = m_Known2geneidsMap.get(tid2);
			if (geneid1 == null || geneid2 == null)
				continue;
			if (geneid1.equals(geneid2))
				continue;
			if (complex_tid2tid_set.contains(tid2tid))
				num_within_complex++;
			if (complex_tid_set.contains(tid1) && complex_tid_set.contains(tid2)){
				num_withComplex++;
			}
		}
		debug_println("# interactions within complex = " + num_within_complex / 2 );
		//debug_println("# all interactions = " + GSP_IntAct_set.size() / 2 );
		debug_println("# all interactions with complex isoforms = " + num_withComplex / 2 );
		debug_println("(# interactions within complex) / (# all interactions) = " 
				+ (double)num_within_complex / num_withComplex);
	}
	
	void ComplexAnalysis_protein() throws Exception {
		Set<String> complex_pid2pid_set = FileUtil.getFileAsSet("../complex_pid2pid_set_" + GSP_one2one_flag);
		Set<String> complex_pid_set = new HashSet<String>();
		for (String s: complex_pid2pid_set){
			String[] items = s.split(":"); 
			complex_pid_set.add(items[0]);
			complex_pid_set.add(items[1]);
		}
		debug_println(" #complex cover proteins  = " + complex_pid_set.size() );
		getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		int num_within_complex = 0;
		int num_withComplex = 0;
		for (String pid2pid: m_pidInt_set){
			String[] items = pid2pid.split(":");
			if (complex_pid2pid_set.contains(pid2pid))
				num_within_complex++;
			if (complex_pid_set.contains(items[0]) && complex_pid_set.contains(items[1])){
				num_withComplex++;
			}
		}
		//complex_pid2pid_set.retainAll(m_pidInt_set);
		debug_println("# interactions within complex = " + num_within_complex / 2 );
		//debug_println(" # all interactions = " + m_pidInt_set.size() / 2 );
		debug_println("# all interactions with complex proteins = " + num_withComplex / 2 );
		debug_println("(# interactions within complex) / (# all interactions) = " 
				+ (double)num_within_complex / num_withComplex);
	}
	
	// TODO: done. 
	void ComplexAnalysis_random() throws Exception { //LR_th = 52794, 64975
		Set<String> complex_tid2tid_set = FileUtil.getFileAsSet("../complex_tid2tid_set_" + GSP_one2one_flag);	
		Set<String> complex_tid_set = new HashSet<String>();
		for (String s: complex_tid2tid_set){
			String[] items = s.split(":"); 
			complex_tid_set.add(items[0]);
			complex_tid_set.add(items[1]);
		}
		debug_println(" # complex cover isoforms  = " + complex_tid_set.size() );
		
		//List<String> isoform_tid2tid_list = getFileIsoformIntList(); 
		//List<String> isoform_tid2tid_list = new ArrayList<String>(); 
		Set<String> tid_set = new HashSet<String>();
		List<String> lines = FileUtil.getFileAsList("out_interaction_features.log");
		for (String s: lines){
			String[] items = s.split("\t"); 
			tid_set.add(items[0]);
			tid_set.add(items[3]);
		}
		debug_println(" # cover isoforms from file = " + tid_set.size() );

		//List<String> random_tid2tid_list = getRandomList(isoform_tid2tid_list, (int)LR_th*2);
		
		Set<String> random_tid_set = new HashSet<String>();
		List<String> tid_list = new ArrayList<String>(tid_set);
		int num_random_Int = 0;
		int num_withinComplex = 0;
		int num_withComplex = 0;
		while (num_withComplex < (int)LR_th) {
			int random_i = random.nextInt(tid_list.size());
			String tid1 = tid_list.get(random_i);
			int random_j = random.nextInt(tid_list.size());
			String tid2 = tid_list.get(random_j);
			String geneid1 = m_Known2geneidsMap.get(tid1);
			String geneid2 = m_Known2geneidsMap.get(tid2);
			if (geneid1 == null || geneid2 == null)
				continue;
			if (tid1.equals(tid2)){
				continue;
			}
			num_random_Int++;
			if (num_random_Int % 1000 == 0)
				debug_println(" # num_random_Int "+ num_random_Int);
			random_tid_set.add(tid1);
			random_tid_set.add(tid2);
			String tid2tid = tid1 + ":" + tid2;
			String reverse_tid2tid = tid2 + ":" + tid1;
			if (complex_tid_set.contains(tid1) && complex_tid_set.contains(tid2)){
				num_withComplex++;
			}
			if (complex_tid2tid_set.contains(tid2tid) || complex_tid_set.contains(reverse_tid2tid)){
				num_withinComplex++;
			}
		}		
		debug_println(" # random interactions = " + num_random_Int );
		debug_println("# random Covered transcripts  = " + random_tid_set.size() );
		debug_println("# interactions within complex = " + num_withinComplex );
		debug_println("# all interactions with complex isoforms = " + num_withComplex);
		debug_println("(# interactions within complex) / (# all interactions) = "
				+ (double)num_withinComplex / num_withComplex);
	}
	
	// TODO: done. 
	void ComplexAnalysis_random_isoform() throws Exception {
		Set<String> complex_tid2tid_set = FileUtil.getFileAsSet("../complex_tid2tid_set_" + GSP_one2one_flag);
		Set<String> complex_tid_set = new HashSet<String>();
		for (String s: complex_tid2tid_set){
			String[] items = s.split(":"); 
			complex_tid_set.add(items[0]);
			complex_tid_set.add(items[1]);
		}
		debug_println(" #complex cover isoforms  = " + complex_tid_set.size() );
		
		//List<String> isoform_tid2tid_list = getAllIsoformIntList(complex_tid_set); 
		getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		List<String> isoform_tid2tid_list = new ArrayList<String>(); 
		Set<String> tid_set = new HashSet<String>();
		for (String tid2tid: m_tidInt2pidInt.keySet()){
			String[] items = tid2tid.split(":"); 
			tid_set.add(items[0]);
			tid_set.add(items[1]);
		}
		debug_println(" # cover isoforms from database = " + tid_set.size() );
		debug_println(" # all interactions within cover isoforms = " + isoform_tid2tid_list.size() ); 
		
		Set<String> random_tid_set = new HashSet<String>();
		List<String> tid_list = new ArrayList<String>(tid_set);
		int num_random_Int = 0;
		int num_withinComplex = 0;
		int num_withComplex = 0;
		while (num_withComplex < (int)LR_th) {
			int random_i = random.nextInt(tid_list.size());
			String tid1 = tid_list.get(random_i);
			int random_j = random.nextInt(tid_list.size());
			String tid2 = tid_list.get(random_j);
			String geneid1 = m_Known2geneidsMap.get(tid1);
			String geneid2 = m_Known2geneidsMap.get(tid2);
			if (geneid1 == null || geneid2 == null)
				continue;
			if (tid1.equals(tid2)){
				continue;
			}
			num_random_Int++;
			if (num_random_Int % 1000 == 0)
				debug_println(" # num_random_Int "+ num_random_Int);
			random_tid_set.add(tid1);
			random_tid_set.add(tid2);
			String tid2tid = tid1 + ":" + tid2;
			String reverse_tid2tid = tid2 + ":" + tid1;
			if (complex_tid_set.contains(tid1) && complex_tid_set.contains(tid2)){
				num_withComplex++;
			}
			if (complex_tid2tid_set.contains(tid2tid) || complex_tid_set.contains(reverse_tid2tid)){
				num_withinComplex++;
			}
		}		
		debug_println(" # random interactions = " + num_random_Int );
		debug_println("# random Covered transcripts  = " + random_tid_set.size() );
		debug_println("# interactions within complex = " + num_withinComplex );
		debug_println("# all interactions with complex isoforms = " + num_withComplex);
		debug_println("(# interactions within complex) / (# all interactions) = "
				+ (double)num_withinComplex / num_withComplex);
	
		/*List<String> random_tid2tid_list = getRandomList(isoform_tid2tid_list, (int)LR_th*2);
		debug_println(" ## random interactions within cover isoforms = " + random_tid2tid_list.size() ); 
		Set<String> random_tid2tid_set = new HashSet<String>(random_tid2tid_list);
		debug_println(" ## random interactions within cover isoforms = " + random_tid2tid_set.size() ); 
		debug_println(" # random interactions = " + random_tid2tid_set.size() / 2 );
		Set<String> random_tid_set = new HashSet<String>();
		for (String s: random_tid2tid_set){
			String[] items = s.split(":"); 
			random_tid_set.add(items[0]);
			random_tid_set.add(items[1]);
		}
		debug_println(" # random cover isoforms  = " + complex_tid_set.size() );
		
		int num_withComplex = 0;
		for (String tid2tid: random_tid2tid_set){
			String[] items = tid2tid.split(":");
			if (complex_tid_set.contains(items[0]) && complex_tid_set.contains(items[1])){
				num_withComplex++;
			}
		}
		complex_tid2tid_set.retainAll(random_tid2tid_set);
		debug_println("# interactions within complex = " + complex_tid2tid_set.size() / 2 );
		debug_println("# all interactions with complex isoforms = " + num_withComplex / 2 );
		debug_println("(# interactions within complex) / (# all interactions) = "
				+ (double)complex_tid2tid_set.size() / num_withComplex);*/
	}

	// TODO: done. 
	void ComplexAnalysis_random_protein() throws Exception {
		Set<String> complex_pid2pid_set = FileUtil.getFileAsSet("../complex_pid2pid_set_" + GSP_one2one_flag);
		Set<String> complex_pid_set = new HashSet<String>();
		for (String s: complex_pid2pid_set){
			String[] items = s.split(":"); 
			complex_pid_set.add(items[0]);
			complex_pid_set.add(items[1]);
		}
		debug_println(" #complex cover proteins  = " + complex_pid_set.size() );
		
		//List<String> protein_pid2pid_list = getAllProteinIntList(complex_pid_set); 
		getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		Set<String> pid_set = new HashSet<String>();
		for (String pid2pid: m_pidInt_set){
			String[] items = pid2pid.split(":"); 
			pid_set.add(items[0]);
			pid_set.add(items[1]);
		}
		debug_println(" # cover proteins from database = " + pid_set.size() );
		debug_println(" # all interactions within cover proteins = " + m_pidInt_set.size() ); 
		
		Set<String> random_pid_set = new HashSet<String>();
		List<String> pid_list = new ArrayList<String>(pid_set);
		int num_random_Int = 0;
		int num_withinComplex = 0;
		int num_withComplex = 0;
		while (num_withComplex < (int)LR_th) {
			int random_i = random.nextInt(pid_list.size());
			String pid1 = pid_list.get(random_i);
			int random_j = random.nextInt(pid_list.size());
			String pid2 = pid_list.get(random_j);
			if (pid1.equals(pid2)){
				continue;
			}
			num_random_Int++;
			if (num_random_Int % 1000 == 0)
				debug_println(" # num_random_Int "+ num_random_Int);
			random_pid_set.add(pid1);
			random_pid_set.add(pid2);
			String pid2pid = pid1 + ":" + pid2;
			String reverse_pid2pid = pid2 + ":" + pid1;
			if (complex_pid_set.contains(pid1) && complex_pid_set.contains(pid2)){
				num_withComplex++;
			}
			if (complex_pid2pid_set.contains(pid2pid) || complex_pid_set.contains(reverse_pid2pid)){
				num_withinComplex++;
			}
		}		
		debug_println(" # random interactions = " + num_random_Int );
		debug_println("# random Covered proteins  = " + random_pid_set.size() );
		debug_println("# interactions within complex = " + num_withinComplex );
		debug_println("# all interactions with complex proteins = " + num_withComplex);
		debug_println("(# interactions within complex) / (# all interactions) = "
				+ (double)num_withinComplex / num_withComplex);
		
		/*List<String> random_pid2pid_list = getRandomList(protein_pid2pid_list, (int)LR_th*2);
		debug_println(" ## random interactions within cover proteins = " + random_pid2pid_list.size() ); 
		Set<String> random_pid2pid_set = new HashSet<String>(random_pid2pid_list);
		debug_println(" ## random interactions within cover proteins = " + random_pid2pid_set.size() ); 
		debug_println(" # random interactions = " + random_pid2pid_set.size() / 2 );
		Set<String> random_pid_set = new HashSet<String>();
		for (String s: random_pid2pid_set){
			String[] items = s.split(":"); 
			random_pid_set.add(items[0]);
			random_pid_set.add(items[1]);
		}
		debug_println(" # random cover proteins  = " + complex_pid_set.size() );

		int num_withComplex = 0;
		for (String pid2pid: random_pid2pid_set){
			String[] items = pid2pid.split(":");
			if (complex_pid_set.contains(items[0]) && complex_pid_set.contains(items[1])){
				num_withComplex++;
			}
		}		
		complex_pid2pid_set.retainAll(random_pid2pid_set);
		debug_println("# interactions within complex = " + complex_pid2pid_set.size() / 2 );
		debug_println("# all interactions with complex proteins = " + num_withComplex / 2 );
		debug_println("(# interactions within complex) / (# all interactions) = "
				+ (double)complex_pid2pid_set.size() / num_withComplex);*/
	}
	
	List<String> getFileIsoformIntList() throws Exception {
		List<String> isoform_tid2tid_list = new ArrayList<String>(); 
		Set<String> tid_set = new HashSet<String>();
		List<String> lines = FileUtil.getFileAsList("out_interaction_features.log");
		for (String s: lines){
			String[] items = s.split("\t"); 
			tid_set.add(items[0]);
			tid_set.add(items[3]);
		}
		debug_println(" # cover isoforms from file = " + tid_set.size() );
		int i = 0;
		for(String tid1: tid_set) {
			i++;
			if (i % 1000 == 0)
				debug_println("#tid "+ i);
			for(String tid2: tid_set) {
				if (tid1.equals(tid2)){
					continue;
				}
				String tid2tid = tid1 + ":" + tid2;
				isoform_tid2tid_list.add(tid2tid);
			}
		}		
		return isoform_tid2tid_list;
	}
	
	List<String> getAllIsoformIntList(Set<String> complex_tid_set) throws Exception {
		getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		List<String> isoform_tid2tid_list = new ArrayList<String>(); 
		Set<String> tid_set = new HashSet<String>();
		for (String tid2tid: m_tidInt2pidInt.keySet()){
			String[] items = tid2tid.split(":"); 
			tid_set.add(items[0]);
			tid_set.add(items[1]);
		}
		debug_println(" # cover isoforms from database = " + tid_set.size() ); 
		int i = 0;
		for(String tid1: tid_set) {
			i++;
			if (i % 1000 == 0)
				debug_println("#tid "+ i);
			for(String tid2: tid_set) {
				if (tid1.equals(tid2)){
					continue;
				}
				if (!complex_tid_set.contains(tid1) 
						&& !complex_tid_set.contains(tid2))
					continue;
				String tid2tid = tid1 + ":" + tid2;
				isoform_tid2tid_list.add(tid2tid);
			}
		}		
		return isoform_tid2tid_list;
	}
	
	List<String> getAllProteinIntList(Set<String> complex_pid_set) throws Exception {
		getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		List<String> protein_pid2pid_list = new ArrayList<String>(); 
		Set<String> pid_set = new HashSet<String>();
		for (String pid2pid: m_pidInt_set){
			String[] items = pid2pid.split(":"); 
			pid_set.add(items[0]);
			pid_set.add(items[1]);
		}
		debug_println(" # cover proteins from database = " + pid_set.size() ); 
		int i = 0;
		for(String pid1: pid_set) {
			i++;
			if (i % 1000 == 0)
				debug_println("#pid "+ i);
			for(String pid2: pid_set) {
				if (pid1.equals(pid2)){
					continue;
				}
				if (!complex_pid_set.contains(pid1) 
						&& !complex_pid_set.contains(pid2))
					continue;
				String pid2pid = pid1 + ":" + pid2;
				protein_pid2pid_list.add(pid2pid);
			}
		}			
		return protein_pid2pid_list;
	}
	
	List<String> getRandomList(List<String> tid2tid_list, int size) {
		List<String> result = new ArrayList<String>();
		while (result.size() < size) {
			int i = random.nextInt(tid2tid_list.size());
			String tid2tid = tid2tid_list.get(i);
			if (result.contains(tid2tid))
				continue;
			String[] items = tid2tid.split(":"); 
			String reverse_tid2tid = items[1] + ":" + items[0];
			result.add(tid2tid);
			result.add(reverse_tid2tid);
		}
		return result;
	}
	
	// TODO: done. 
	void ComplexAnalysis(List<Double> CoefList) throws Exception { 
		for (double c: CoefList){
			System.out.println(c);
		}
		getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", GSP_one2one_flag);
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
		BufferedWriter bw1 = new BufferedWriter(new FileWriter("out_interaction_features.log"));
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
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0723_2009");
			Tid2goname = getTid2gonameUniProt("uniprot2go_0405_2011_BP");
		}	
		Map<String, Double> HomologyRatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		Set<String> complex_tid2tid_set = FileUtil.getFileAsSet("../complex_tid2tid_set_" + GSP_one2one_flag);
		Set<String> complex_tid_set = new HashSet<String>();
		for (String s: complex_tid2tid_set){
			String[] items = s.split(":"); 
			complex_tid_set.add(items[0]);
			complex_tid_set.add(items[1]);
		}
		debug_println(" # complex_tid_set = " + complex_tid_set.size() );
		int num_all_interaction = 0;
		int num_within_complex = 0;
		int num_withComplex = 0;
		for(int i=0; i < tid_list.size(); i++) {
			if (i % 1000 == 0)
				debug_println("#tid "+ i);
			String tid1 = tid_list.get(i);
			for(int j=0; j < tid_list.size(); j++) {
				String tid2 = tid_list.get(j);
				if (tid1.equals(tid2)){
					continue;
				}
				String geneid1 = m_Known2geneidsMap.get(tid1);
				String geneid2 = m_Known2geneidsMap.get(tid2);
				if (geneid1.equals(geneid2))
					continue;
				String tid2tid = tid1 + ":" + tid2;
				/*if (!GSP_IntAct_set.contains(tid2tid)
						&& !((GSN_set.tids_membrane.contains(tid1) && 
							  GSN_set.tids_nucleus.contains(tid2)) ||
							 (GSN_set.tids_membrane.contains(tid2) && 
							  GSN_set.tids_nucleus.contains(tid1)))
							  ){
					continue;					
				}*/
				
				String R_str = getExpScore_str(i, j, Profile_list);
				String[] items = R_str.split(",");
				//d1, d2, d3, d4 ,d5, d7, d9, d10, d11
				double exp_d1_score = Double.parseDouble(items[0]);
				double exp_d2_score = Double.parseDouble(items[1]);
				double exp_d3_score = Double.parseDouble(items[2]);
				double exp_d4_score = Double.parseDouble(items[3]);
				double exp_d5_score = Double.parseDouble(items[4]);
				double exp_d7_score = Double.parseDouble(items[5]);
				double exp_d9_score = Double.parseDouble(items[6]);
				double exp_d10_score = Double.parseDouble(items[7]);
				double exp_d11_score = Double.parseDouble(items[8]);
				double domain_score = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double GO_score = getGoScore(tid1, tid2, Tid2goname);
				double ortholog_score = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				double LR = 
					exp_d1_score * CoefList.get(0) + exp_d2_score * CoefList.get(1) +
					exp_d3_score * CoefList.get(2) + exp_d4_score * CoefList.get(3) +
				    exp_d5_score * CoefList.get(4) + exp_d7_score * CoefList.get(5) + 
				    exp_d9_score * CoefList.get(6) + exp_d10_score * CoefList.get(7) + 
				    exp_d11_score * CoefList.get(8) + 
				    domain_score * CoefList.get(9) + 
				    GO_score * CoefList.get(10) + 
				    ortholog_score * CoefList.get(11);	
				
				//Calculate (# interactions within complex) / (# all interactions)
				if (LR > LR_th){
					outInteractionFile2(LR, tid1, tid2, geneid1, geneid2, bw1);
					num_all_interaction++;
					
					if (complex_tid2tid_set.contains(tid2tid) 
							&& m_tidInt2pidInt.keySet().contains(tid2tid)
							)
						num_within_complex++;
					if (complex_tid_set.contains(tid1) && complex_tid_set.contains(tid2)
							&& m_all_knowns_set.contains(tid1) && m_all_knowns_set.contains(tid2)
							)
						num_withComplex++;
				}
			}
		}
		//# interactions within complex = 5952 
		//# all interactions = 507573
		//# all interactions with complex proteins = 
		//(# interactions within complex) / (# all interactions) = 
		debug_println(" # interactions within complex = " + num_within_complex / 2 );
		debug_println(" # all interactions = " + num_all_interaction / 2 );
		debug_println(" # all interactions with complex proteins = " + num_withComplex / 2 );
		debug_println(" (# interactions within complex) / (# all interactions) = " 
				+ (double)num_within_complex / num_withComplex);
		bw1.close();
	}
	
	static void ExcuteR() throws Exception {
		String cmdline = "nohup ./R.sh &";
		Process ps = Runtime.getRuntime().exec(cmdline);
		try {
			ps.waitFor();
		} catch (InterruptedException x) {
		}
	}
	
	static List<Double> getCoefList() throws Exception{
		List<Double> coef_list = new ArrayList<Double>();
		BufferedReader br = new BufferedReader(new FileReader("./temp/coef.out"));
		String line;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("x")||(line.startsWith("(Intercept)")))
				continue;
			String[] items = line.split(" ");
			double coef = Double.parseDouble(items[1]);
			coef_list.add(coef);
		}
		br.close();
		debug_println("getCoefList() CoefList "+ coef_list.size());
		return coef_list;	
	}
	
	void prepare_all_arr() throws Exception{
		//Arr all_arr = new Arr();
		
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
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("Regression_IsoIntAct.csv"));
		bw.write("Int_score,exp_d1_score,exp_d2_score,exp_d3_score,exp_d4_score,exp_d5_score,exp_d7_score," +
				"exp_d9_score,exp_d10_score,exp_d11_score,exp_SRP002079,domain_score,GO_score,ortholog_score");
		bw.write("\n");
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
					 && !Train_GSN_set.contains(tid2tid)
							  ){
					continue;					
				}
			
				String R_str = getExpScore_str(i, j, Profile_list);
				double maxRatioD = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double SSBP = getGoScore(tid1, tid2, Tid2goname);
				double Class = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				if (Train_GSP_set.contains(tid2tid)){		
					bw.write("1," + R_str +  maxRatioD + "," + SSBP + "," + Class); 
					bw.write("\n");
					num_pos++;
				}
				if (Train_GSN_set.contains(tid2tid)){
					//if (num_pos > num_neg) {
						bw.write("0," + R_str +  maxRatioD + "," + SSBP + "," + Class); 
						bw.write("\n");
						num_neg++;
					//}
				}
			}
		}
		bw.close();
		debug_println("  num_pos " + num_pos);
		debug_println("  num_neg " + num_neg);
	}
	
	/**
	* 12 features: Use 9 RNA-seq + domain + GO + ortholog. (d1, d2, d3, d4 ,d5, d7, d9, d10, d11)
	* 11 features: Use 9 RNA-seq + domain + ortholog. (d1, d2, d3, d4 ,d5, d7, d9, d10, d11)
	* 8 features: Use 5 RNA-seq + domain + GO + ortholog. (d1, d2, d7, d9, d11)
	* 7 features: Use 5 RNA-seq + domain + ortholog. (d1, d2, d7, d9, d11)
	*/
	
	//uniprot2go_0405_2011_BP: Biological process GO for GO score.
	//uniprot2go_0405_2011: All GO for GSN construction.
	
	void validation(List<Double> CoefList) throws Exception { 
		for (double c: CoefList){
			System.out.println(c);
		}

		debug_println("Known2geneid " + m_Known2geneidsMap.size());
		BufferedWriter bw1 = new BufferedWriter(new FileWriter("interaction.log"));
		int gsp = 0;
		int gsn = 0;
		int test_pos = 0;
		int test_neg = 0;
		int total = 0;
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
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0723_2009");
			Tid2goname = getTid2gonameUniProt("uniprot2go_0405_2011_BP");
		}	
		Map<String, Double> HomologyRatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);

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
				
				//String R_str = getExpScore_str(i, j, Profile_list);
				//String[] items = R_str.split(",");
				double R_ary[] = getExpScore_ary(i, j, Profile_list);
				//d1, d2, d3, d4 ,d5, d7, d9, d10, d11
				double exp_d1_score = R_ary[0];
				double exp_d2_score = R_ary[1];
				double exp_d3_score = R_ary[2];
				double exp_d4_score = R_ary[3];
				double exp_d5_score = R_ary[4];
				double exp_d7_score = R_ary[5];
				double exp_d9_score = R_ary[6];
				double exp_d10_score = R_ary[7];
				double exp_d11_score = R_ary[8];
				double domain_score = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double GO_score = getGoScore(tid1, tid2, Tid2goname);
				double ortholog_score = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				double y_score = 
					exp_d1_score * CoefList.get(0) + exp_d2_score * CoefList.get(1) +
					exp_d3_score * CoefList.get(2) + exp_d4_score * CoefList.get(3) +
				    exp_d5_score * CoefList.get(4) + exp_d7_score * CoefList.get(5) + 
				    exp_d9_score * CoefList.get(6) + exp_d10_score * CoefList.get(7) + 
				    exp_d11_score * CoefList.get(8) + 
				    domain_score * CoefList.get(9) + 
				    GO_score * CoefList.get(10) + 
				    ortholog_score * CoefList.get(11);
				
				/*double z = exp_d2_score * CoefList.get(0) + exp_d3_score * CoefList.get(1) 
						+ exp_d4_score * CoefList.get(2) + exp_d5_score * CoefList.get(3) 
						+ exp_d9_score * CoefList.get(4) + domain_score * CoefList.get(5)
						+ GO_score * CoefList.get(6) + ortholog_score * CoefList.get(7);
				double y_score = 1 / (1 + Math.exp(-z));*/
	
				outInteractionFile(y_score, tid1, tid2, bw1);
				total++;
				if (Train_GSP_set.contains(tid2tid)){
					gsp++;					
				}
				if (Train_GSN_set.contains(tid2tid)){
					gsn++;
				}
				if (Test_GSP_set.contains(tid2tid)){
					test_pos++;
				}
				if (Test_GSN_set.contains(tid2tid)){
					test_neg++;
				}
			}
		}
		debug_println("  #gsp " + gsp);
		debug_println("  #gsn " + gsn);		
		debug_println("  #test_pos " + test_pos);
		debug_println("  #test_neg " + test_neg);
		debug_println("  #total " + total);
		bw1.close();
	}
	
	void out_validation(List<Double> CoefList) throws Exception { 
		for (double c: CoefList){
			System.out.println(c);
		}
		Set<String> GSP_iRef_set = getGSPsetiRef("/home/jimliu/IntAct/iRefWeb_lite.txt");
		debug_println("GSP_iRef_set " + GSP_iRef_set.size());
		Set<String> GSP_IntAct_set = getGSPsetIntActIso("/home/ting/IntAct/2011_03_09/psimitab/intact.txt", false);
		debug_println("Known2geneid " + m_Known2geneidsMap.size());
		//Map<String, String> m_Known2protein = GeneAnnot.getKnown2protein(m_conn);
		//Map<String, String> m_Known2protein = GeneAnnot.getKnown2protein(hg18_to_uniprot);
		Map<String, String> m_Known2protein = null;
		if (hg18_to_uniprot == null){
			m_Known2protein = GeneAnnot.getKnown2protein(m_conn);
		}else{
			m_Known2protein = GeneAnnot.getKnown2protein(hg18_to_uniprot);
		}
		BufferedWriter bw1 = new BufferedWriter(new FileWriter("out_interaction_features.log"));
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
			//Tid2goname = getTid2gonameUniProt("uniprot2go_0723_2009");
			Tid2goname = getTid2gonameUniProt("uniprot2go_0405_2011_BP");
		}	
		Map<String, Double> HomologyRatioMap = getHomologyRatioMap(HomologyRatioMap_fileName);

		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		
		bw1.write("Tid1\tGeneId1\tGeneSym1\tProteinId1\t" +
		          "Tid2\tGeneId2\tGeneSym2\tProteinId2\t" +
		          "Score\tIntAct\tiRef\texp\tdomain\tGO\tortholog\t" +
				  "exp_d1\texp_d2\texp_d3\texp_d4\texp_d5\t" +
				  "exp_d7\texp_d9\texp_d10\texp_d11\texp_SRP002079\n");
		for(int i=0; i < tid_list.size(); i++) {
			if (i % 100 == 0)
				debug_println("#tid "+ i);
			String tid1 = tid_list.get(i);
			for(int j=0; j < tid_list.size(); j++) {
				String tid2 = tid_list.get(j);
				if (tid1.equals(tid2)){
					continue;
				}
				String geneid1 = m_Known2geneidsMap.get(tid1);
				String geneid2 = m_Known2geneidsMap.get(tid2);
				if ((geneid1 != null) && (geneid1 != null)) {
					if (geneid1.equals(geneid2))
						continue;
				}
				//String tid2tid = tid1 + ":" + tid2;
				/*if (!GSP_IntAct_set.contains(tid2tid)
						&& !((GSN_set.tids_membrane.contains(tid1) && 
							  GSN_set.tids_nucleus.contains(tid2)) ||
							 (GSN_set.tids_membrane.contains(tid2) && 
							  GSN_set.tids_nucleus.contains(tid1)))
							  ){
					continue;					
				}*/
				/*if (!Train_GSP_set.contains(tid2tid) 
						&& !Train_GSN_set.contains(tid2tid) 
						&& !Test_GSP_set.contains(tid2tid) 
						&& !Test_GSN_set.contains(tid2tid) 
							  ){
					continue;					
				}*/
				
				//String R_str = getExpScore_str(i, j, Profile_list);
				//String[] items = R_str.split(",");
				double R_ary[] = getExpScore_ary(i, j, Profile_list);
				//d1, d2, d3, d4 ,d5, d7, d9, d10, d11, SRP002079
				double exp_d1_score = R_ary[0];
				double exp_d2_score = R_ary[1];
				double exp_d3_score = R_ary[2];
				double exp_d4_score = R_ary[3];
				double exp_d5_score = R_ary[4];
				double exp_d7_score = R_ary[5];
				double exp_d9_score = R_ary[6];
				double exp_d10_score = R_ary[7];
				double exp_d11_score = R_ary[8];
				double exp_SRP002079_score = R_ary[9];
				
				double exp_score = MathUtil.getAver(R_ary);
				double domain_score = getDomainScore(tid1, tid2, tid2Pfam, RatioMap);
				double GO_score = getGoScore(tid1, tid2, Tid2goname);
				double ortholog_score = getOrthologScore(tid1, tid2, HomologyRatioMap);
				
				double y_score = 
					exp_d1_score * CoefList.get(0) + exp_d2_score * CoefList.get(1) +
					exp_d3_score * CoefList.get(2) + exp_d4_score * CoefList.get(3) +
				    exp_d5_score * CoefList.get(4) + exp_d7_score * CoefList.get(5) + 
				    exp_d9_score * CoefList.get(6) + exp_d10_score * CoefList.get(7) + 
				    exp_d11_score * CoefList.get(8) + exp_SRP002079_score * CoefList.get(9) + 
				    domain_score * CoefList.get(10) + 
				    GO_score * CoefList.get(11) + 
				    ortholog_score * CoefList.get(12);
				
				//outInteractionFile(y_score, tid1, tid2, bw1);
				
				//Format: tid1 \t gene_id1 \t gene_sym1 \t protein_id1 \t 
				//        tid2 \t gene_id2 \t gene_sym2 \t protein_id2 \t 
				//        score \t IntAct interaction \t iRef interaction
				//        exp1 \t exp2 \t ... \t domain \t ...
				if (y_score > LR_th){
					
					bw1.write(tid1 + "\t ");
					if (geneid1 != null) {
						bw1.write(geneid1 + "\t " + GeneAnnot.getGeneSym(geneid1) + "\t");
					}else{
						bw1.write("null\tnull\t");
					}
					String proteinid1 = m_Known2protein.get(tid1);
					if (proteinid1 != null) {
						bw1.write(proteinid1 + "\t ");
					}else{
						bw1.write("null\t");
					}
					
					bw1.write(tid2 + "\t ");
					if (geneid2 != null) {
						bw1.write(geneid2 + "\t " + GeneAnnot.getGeneSym(geneid2) + "\t");
					}else{
						bw1.write("null\tnull\t");
					}
					String proteinid2 = m_Known2protein.get(tid2);
					if (proteinid2 != null) {
						bw1.write(proteinid2 + "\t ");
					}else{
						bw1.write("null\t");
					}
					
					if (m_tidInt2pidInt.get(tid1 + ":" + tid2) != null){
						bw1.write(y_score + "\t" +m_tidInt2pidInt.get(tid1 + ":" + tid2)+ "\t");
					}else{
						bw1.write(y_score + "\t" +"null"+ "\t");
					}
					
					if (proteinid1 != null && proteinid1.contains("-")) {
						proteinid1 = proteinid1.substring(0, proteinid1.indexOf('-'));
					}
					if (proteinid2 != null && proteinid2.contains("-")) {
						proteinid2 = proteinid2.substring(0, proteinid2.indexOf('-'));
					}
					if (GSP_iRef_set.contains(proteinid1 + ":" + proteinid2)){
						bw1.write(proteinid1 + ":" + proteinid2 + "\t");
					}else{
						bw1.write("null"+ "\t");
					}
					
					bw1.write(exp_score + "\t" + domain_score + "\t" + GO_score + "\t" + ortholog_score + "\t"
							+ exp_d1_score + "\t" + exp_d2_score + "\t" + exp_d3_score + "\t" 
							+ exp_d4_score + "\t" + exp_d5_score + "\t" + exp_d7_score + "\t" 
							+ exp_d9_score + "\t" + exp_d10_score + "\t" + exp_d11_score + "\t" 
							+ exp_SRP002079_score + "\n");
				}
			}
		}
		bw1.close();
	}
	
	static Set<String> getGSPsetiRef(String fileName) throws Exception {
		debug_println("getGSPsetiRef()");
		List<String> lines = FileUtil.getFileAsList(fileName);
		Set<String> m_pidInt_set2 = new HashSet<String>();
		int i = 0;
		for (String line: lines) {
			if (line.startsWith("interaction_id")) {
				continue;
			}
			String[] items = line.split("\t");
			String protein1 = items[1];
			String protein2 = items[2];
			protein1 = protein1.substring(0, 6);
			protein2 = protein2.substring(0, 6);
			m_pidInt_set2.add(protein1 +":" + protein2);
			m_pidInt_set2.add(protein2 +":" + protein1);
			i++;
			if (i<=50){
				debug_println(i + " " + protein1 +":" + protein2);
			}
		}
		return m_pidInt_set2;
	}

	static void outInteractionFile(double LR, String tid1, String tid2, BufferedWriter bw1) throws Exception{
		//Format: tid1 \t gene_id1 \t tid2 \t gene_id2 \t score
		if (LR > LR_th){
			bw1.write(tid1 + "\t");
			String geneid1 = m_Known2geneidsMap.get(tid1);
			if (geneid1 != null) {
				bw1.write(geneid1 + "\t");
			}else{
				bw1.write("null\t");
			}
			String geneid2 = m_Known2geneidsMap.get(tid2);
			
			bw1.write(tid2 + "\t");
			if (geneid2 != null) {
				bw1.write(geneid2 + "\t");
			}else{
				bw1.write("null\t");
			}
			bw1.write(LR + "\n");
		}
	}
	
	static void outInteractionFile2(double LR, String tid1, String tid2, String geneid1, String geneid2, BufferedWriter bw1) 
	//Use the threshold 1.0394, output the III prediction.
	//Format: tid1 \t gene_id1 \t gene_sym1 \t tid2 \t gene_id2 \t gene_sym2 \t score \t IntAct interaction
	throws Exception{
		bw1.write(tid1 + "\t");
		if (geneid1 != null) {
			bw1.write(geneid1 + "\t" + GeneAnnot.getGeneSym(geneid1) + "\t");
		}else{
			bw1.write("null\tnull\t");
		}
		
		bw1.write(tid2 + "\t");
		if (geneid2 != null) {
			bw1.write(geneid2 + "\t" + GeneAnnot.getGeneSym(geneid2) + "\t");
		}else{
			bw1.write("null\tnull\t");
		}
		
		if (m_tidInt2pidInt.get(tid1 + ":" + tid2) != null){
			bw1.write(LR + "\t" +m_tidInt2pidInt.get(tid1 + ":" + tid2)+ "\n");
		}else{
			bw1.write(LR + "\t" +"null"+ "\n");
		}
	}
	
	static String getExpScore_str(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		//double R = 0;
		String R_str = "";
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			double R = 0;
			if (exp[i] != null && exp[j] != null) {
				//R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
				R = MathUtil.getPearson(exp[i], exp[j]);
			}
			R_str = R_str + R + ",";
		}
		return R_str;
	}
	
	static double[] getExpScore_ary(int i, int j, List<ExpProfile> profile_list) 
	throws Exception {
		//double R = 0;
		double[] R_ary = new double[10];
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			double R = 0;
			if (exp[i] != null && exp[j] != null) {
				//R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
				R = MathUtil.getPearson(exp[i], exp[j]);
			}
			R_ary[net] = R;
		}
		return R_ary;
	}
	
	/*String getExpScore_n_str(int i, int j, List<ExpProfile> profile_list, Features features) 
	throws Exception {
		double R = 0;
		double n_R = 0;
		String n_R_str = "";
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
			}
			if (net==0) n_R = (R - features.exp_d2_Aver) / features.exp_d2_SD;
			if (net==1) n_R = (R - features.exp_d3_Aver) / features.exp_d3_SD;
			if (net==2) n_R = (R - features.exp_d4_Aver) / features.exp_d4_SD;
			if (net==3) n_R = (R - features.exp_d5_Aver) / features.exp_d5_SD;
			if (net==4) n_R = (R - features.exp_d9_Aver) / features.exp_d9_SD;
			n_R_str = n_R_str + n_R + ",";
		}
		return n_R_str;
	}*/
	
	/*String getExpScore_str(int i, int j, List<ExpProfile> profile_list, Features features) 
	throws Exception {
		double R = 0;
		String R_str = "";
		for (int net = 0; net < profile_list.size(); net++) {
			double[][] exp = profile_list.get(net).exp;
			if (exp[i] != null && exp[j] != null) {
				R = Math.abs(MathUtil.getPearson(exp[i], exp[j]));
			}
			if (net==0) features.exp_d2_list.add(R);
			if (net==1) features.exp_d3_list.add(R);
			if (net==2) features.exp_d4_list.add(R);
			if (net==3) features.exp_d5_list.add(R);
			if (net==4) features.exp_d9_list.add(R);
			R_str = R_str + R + ",";
		}
		return R_str;
	}*/
	

	//TODO: Output complex data with format: proteinID \t complex name
	void outputComplexData() throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter("pid2complex"));
		String complex_fileName = "allComplexes_CORUM.csv";
		List<String> lines = FileUtil.getFileAsList(complex_fileName);
		for (String line: lines) {
			if (line.startsWith("Complex"))
				continue;
			String[] items = line.split(";");
			String complex_name = items[1];
			String organism = items[3];
			if (!organism.equals("Human"))
				continue;
			String subunits = items[4];
			subunits = subunits.replace("(", ""); //P29279,(P60709,P63261)
			subunits = subunits.replace(")", "");
			String[] proteins = subunits.split(","); 
			for (int i = 0; i < proteins.length; i++) {
				bw.write(proteins[i] + "\t" + complex_name + "\n");
			}	
		}
		bw.close();
	}
	
	void ProcessComplexData() throws Exception {
		String complex_fileName = "allComplexes_CORUM.csv";
		List<String> lines = FileUtil.getFileAsList(complex_fileName);
		Map<String,Set<String>> complex2subunits = new HashMap<String, Set<String>>();
		Set<String> complex_pid2pid = new HashSet<String>();
		Set<String> complex_tid2tid = new HashSet<String>();
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		for (String line: lines) {
			if (line.startsWith("Complex"))
				continue;
			String[] items = line.split(";");
			String complex_name = items[1];
			String organism = items[3];
			if (!organism.equals("Human"))
				continue;
			String subunits = items[4];
			subunits = subunits.replace("(", ""); //P29279,(P60709,P63261)
			subunits = subunits.replace(")", "");
			String[] proteins = subunits.split(","); 
			Set<String> subunits_set = new HashSet<String>();
			for (int i = 0; i < proteins.length; i++){
				subunits_set.add(proteins[i]);
			}
			complex2subunits.put(complex_name, subunits_set);
			
			for (String protein1: subunits_set){
				for (String protein2: subunits_set){
					if (protein1.equals(protein2)) {
						continue;
					}

					Set<String> knowns1 = protein2knowns.get(protein1);
					Set<String> knowns2 = protein2knowns.get(protein2);
					if (knowns1 == null || knowns2 == null) {
						continue;
					}
								
					complex_pid2pid.add(protein1+":"+protein2);
					complex_pid2pid.add(protein2+":"+protein1);
					
					if (GSP_one2one_flag) { // check one-to-one
						if (knowns1.size() == 1 && knowns2.size() == 1) { 
							for (String g1 : knowns1) {
								for (String g2 : knowns2) {
									complex_tid2tid.add(g1 + ":" + g2);
									complex_tid2tid.add(g2 + ":" + g1);
								}
							}
						}
					} else {
						if (knowns1.size() >= 1 && knowns2.size() >= 1) {
							for (String g1 : knowns1) {
								for (String g2 : knowns2) {
									complex_tid2tid.add(g1 + ":" + g2);
									complex_tid2tid.add(g2 + ":" + g1);
								}
							}
						}
					}
					
				}
			}
			
		}
		FileUtil.outputMapAsFile(complex2subunits, "complex2subunits_map_" + GSP_one2one_flag);
		FileUtil.outputCollectionAsFile(complex_pid2pid, "complex_pid2pid_set_" + GSP_one2one_flag);
		FileUtil.outputCollectionAsFile(complex_tid2tid, "complex_tid2tid_set_" + GSP_one2one_flag);
	}
	
	void ProcessComplexData2() throws Exception {
		String complex_fileName = "allComplexes_CORUM.csv";
		List<String> lines = FileUtil.getFileAsList(complex_fileName);
		BufferedWriter bw = new BufferedWriter(new FileWriter("complex2subunits_map"));
		for (String line: lines) {
			/*if (line.startsWith("Complex"))
				continue;*/
			String[] items = line.split(";");
			for (int i = 0; i < items.length; i++){
				bw.write(items[i] + "\t");
			}
			bw.write("\n");
		}
		bw.close();
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
