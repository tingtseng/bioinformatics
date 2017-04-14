package edu.nchu.syslab.ppi;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.Set;

import edu.nchu.syslab.rnaseq.BuildExpMatrix;
import edu.nchu.syslab.rnaseq.ExpProfile;
import edu.usc.zhoulab.common.MathUtil;
import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.GeneAnnot;

public class RNAseqAnalysis {
	static ResourceBundle bundler = ResourceBundle.getBundle("common");
	static String DB_URL = bundler.getString("db.url");
	static String DB_USER = bundler.getString("db.user");
	static String DB_PASS = bundler.getString("db.passwd");
	static Connection m_conn = null;
	static String folderName = null;
	static String Tid_fileName = null;
	static boolean m_debug_flag = true;

	static Map<String, String> m_Known2geneidsMap = null;
	static Map<String, Set<String>> m_Geneid2knownsMap = null;
	
	public static void main(String[] args) 
	throws Exception
	{
		if (args.length != 1) {
			System.out.println("Usage: IsoIntAct Tid_fileName");
			// sel_tids.log
			return;
		}
		Tid_fileName = args[0];
		System.out.println("Tid_fileName " + Tid_fileName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new RNAseqAnalysis();
	}

	public RNAseqAnalysis() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		m_Known2geneidsMap = GeneAnnot.getKnown2geneid(m_conn);
		m_Geneid2knownsMap = GeneAnnot.getGeneid2knowns(m_conn);
		List<String> tid_list = FileUtil.getFileAsList(Tid_fileName);
		List<ExpProfile> Profile_list = prepareExpProfile_list(tid_list);
		
		outExpAnalysis(tid_list, Profile_list);
		outNetAnalysis(tid_list, Profile_list);
	}
	
	void outExpAnalysis(List<String> tid_list, List<ExpProfile> profile_list) 
			throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter(Tid_fileName + "_ExpAnalysis.out"));
		for (int i = 0; i < tid_list.size(); i++) {
			String tid = tid_list.get(i);
			double avg_exp = 0.0;
			double sd = 0.0;
			for (int net = 0; net < profile_list.size(); net++) {
				double[][] exp = profile_list.get(net).exp;
				if (exp[i] != null) {
					avg_exp = MathUtil.getAver(exp[i]);
					sd = MathUtil.getSD(exp[i]);
				}
				int sample = net+1;
				if (net>=7) sample+=1;
				String geneid = m_Known2geneidsMap.get(tid);
				String geneSym = GeneAnnot.getGeneSym(geneid);
				bw.write("d" + sample + "\t" + geneSym + "(" + tid + ")\t" + avg_exp + "\t" + sd + "\n");
			}
		}
		bw.close();	
	}
	
	void outNetAnalysis(List<String> tid_list, List<ExpProfile> profile_list)
			throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter(Tid_fileName + "_NetAnalysis.out"));
		for (int i = 0; i < tid_list.size(); i++) {
			String tid1 = tid_list.get(i);
			for (int j = i+1; j < tid_list.size(); j++) {
				String tid2 = tid_list.get(j);
				if (tid1.equals(tid2)) 
					continue;
				double R = 0.0;
				for (int net = 0; net < profile_list.size(); net++) {
					double[][] exp = profile_list.get(net).exp;
					if (exp[i] != null && exp[j] != null) {
						R = MathUtil.getPearson(exp[i], exp[j]);
					}
					int sample = net+1;
					if (net>=7) sample+=1;
					String geneid1 = m_Known2geneidsMap.get(tid1);
					String geneSym1 = GeneAnnot.getGeneSym(geneid1);
					String geneid2 = m_Known2geneidsMap.get(tid2);
					String geneSym2 = GeneAnnot.getGeneSym(geneid2);
					bw.write("d" + sample + "\t" + geneSym1 + "("  + tid1 + ")\t" + geneSym2 + "(" + tid2 + ")\t" + R + "\n");
				}
			}
		}
		bw.close();
	}	

	static List<ExpProfile> prepareExpProfile_list(List<String> tid_list) throws Exception {
		List<ExpProfile> Profile_list = new ArrayList<ExpProfile>();
		for (int i = 1; i <= 11; i++){
			if (i==8 || i==6) continue;
			folderName = "/home/jimliu/RNA-seq/Results/d"+i+"_hg18_known";
			new BuildExpMatrix(folderName, m_conn);
			List<String> sampleList = BuildExpMatrix.getSampleList(folderName);
			Map<String, Map<String, Double>> sample2expMap = BuildExpMatrix.getTransMap(sampleList);
			ExpProfile transExpProfile = BuildExpMatrix.getExpProfile(sampleList, sample2expMap);
			debug_println("transExpProfile d"+i+": " + transExpProfile);
			//transExpProfile.removeLowSD(sd_th); // sd_th 0.1			
			transExpProfile.applyGeneList(tid_list);
			Profile_list.add(transExpProfile);
		}
		debug_println("prepareExpProfile_list() " + Profile_list.size());
		return Profile_list;
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
