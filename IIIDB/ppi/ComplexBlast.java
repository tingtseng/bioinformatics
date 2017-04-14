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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.Set;

import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.GeneAnnot;

public class ComplexBlast {
	static ResourceBundle bundler = ResourceBundle.getBundle("common");
	static String DB_URL = bundler.getString("db.url");
	static String DB_USER = bundler.getString("db.user");
	static String DB_PASS = bundler.getString("db.passwd");
	static Connection m_conn = null;
	
	static String PDB_fileName = "pdb_seqres.txt";
	static boolean m_debug_flag = true;
	
	public static void main(String[] args) throws Exception {
		/*if (args.length != 1) {
			System.out.println("Usage: IsoIntAct PDB_fileName");
			// /home/jimliu/RNA-seq/homology/knownInt_1e-6_0.1.txt
			return;
		}
		PDB_fileName = args[0];*/
		System.out.println("PDB_fileName " + PDB_fileName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
		new ComplexBlast();
	}
	
	public ComplexBlast() throws Exception {
		//parser();
		//StructureAnalysis();
		BindingSiteAnalysis();
	}

	void BindingSiteAnalysis() throws Exception{
		BufferedReader br = new BufferedReader(new FileReader("D:/190/IntAct_webserver/case study document/program/score_2.csv"));
		Set<String> s1_set = new HashSet<String>();
		String line = null;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("tid1")) {
				continue;
			}
			String[] items = line.split(",");
			String pid1 = items[3];
			String pid2 = items[7];
			s1_set.add(pid1+":"+pid2);
			s1_set.add(pid2+":"+pid1);
		}
		br.close();
		
		BufferedReader br1 = new BufferedReader(new FileReader("D:/190/IntAct_webserver/case study document/program/s2.csv"));
		Set<String> s2_set = new HashSet<String>();
		String line2 = null;
		while ((line2 = br1.readLine()) != null) {
			if (line2.startsWith("ProteinA")) {
				continue;
			}
			String[] items = line2.split(",");
			String pid1 = items[0];
			String pid2 = items[1];
			s2_set.add(pid1+":"+pid2);
			s2_set.add(pid2+":"+pid1);
		}
		br1.close();
		
		s1_set.retainAll(s2_set);
		
		for (String s: s1_set){
			System.out.println(s);
		}
	}
	
	
	void parser() throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(PDB_fileName));
		/**
		 * >1a02_N mol:protein length:301  NUCLEAR FACTOR OF ACTIVATED T CELLS
		 * MRGSHHHHHHTDPHASSVPLEWPLSSQSGSYELRI
		 * >1a02_F mol:protein length:56  AP-1 FRAGMENT FOS
		 * MKRRIRRERNKMAAAKSRNRRRELTDTLQAETDQLEDEKSALQTEIANLLKEKEKL
		 * >1a02_J mol:protein length:56  AP-1 FRAGMENT JUN
		 * MKAERKRMRNRIAASKSRKRKLERIARLEEKVKTLKAQNSELASTANMLREQVAQL*/
		String line = null;
		String pdb_id = null;
		boolean protein = false;
		Set<String> seqs;
		Map<String, Set<String>> pdbid2seq = new HashMap<String, Set<String>>();
		while ((line = br.readLine()) != null) {	
			if (line.startsWith(">")) {
				if (line.contains("protein")){
					pdb_id = line.substring(1, line.indexOf("_"));
					protein = true;
				}else{
					protein = false;
				}
			}else{
				if (!protein)
					continue;
				if (pdbid2seq.get(pdb_id) == null){
					seqs = new HashSet<String>();
					seqs.add(line);
					pdbid2seq.put(pdb_id, seqs);
				}else{
					pdbid2seq.get(pdb_id).add(line);
				}		
			}
		}
		br.close();
		
		System.out.println("pdbid2seq.keySet().size() " + pdbid2seq.keySet().size());
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("pdb.fasta"));
		for (String pdbid: pdbid2seq.keySet()) {
			if (pdbid2seq.get(pdbid).size() < 2)
				continue;
			int i = 0;
			for (String seq: pdbid2seq.get(pdbid)){
				i++;
				bw.write(">" + pdbid + "_" + i + "\n");
				bw.write(seq + "\n");
			}
		}
		bw.close();
	}
	
	public static void outputMapAsFile(Map<String, Set<String>> map, String fileName) 
	throws IOException
	{
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for (Object key: map.keySet())
		{
			bw.write(key.toString()+"\t"+map.get(key));
			bw.newLine();
		}
		bw.close();
	}
	
	void StructureAnalysis() throws Exception {
		
		/**
		 * 5gch_1	uc002pvs.1	1	94	172	264	5e-11	63.5	3	36	97
		 * 5gch_1	uc002cro.1	7	91	180	275	8e-11	62.8	3	36	96
		 * 5gch_3	uc002fds.1	1	131	34	164	1e-60	 228	0	109	131
		 * 5gch_3	uc002fdr.1	1	131	34	164	1e-60	 228	0	109	131
		 * 5gch_3	uc002euw.1	1	127	34	161	8e-39	 155	1	72	128 */
		List<String> lines1 = FileUtil.getFileAsList("blastp.out");
		Map<String, Set<String>> tid2pdbid = new HashMap<String, Set<String>>();
		Map<String, Set<String>> tid2pdbid_tmp = new HashMap<String, Set<String>>();
		int j = 0;
		for (String line1: lines1) {
			j++;
			/*if (j % 1000 == 0)
				System.out.println(j);*/
			String[] items = line1.split("\t");
			String tid = items[1];
			String pdbid = items[0];
			Set<String> pdbids = tid2pdbid.get(tid);
			Set<String> pdbids_tmp = tid2pdbid.get(tid);
			if (pdbids == null) {
				pdbids = new HashSet<String>();
				pdbids_tmp = new HashSet<String>();
				pdbids.add(pdbid);
				pdbids_tmp.add(pdbid);
				tid2pdbid.put(tid, pdbids);
				tid2pdbid_tmp.put(tid, pdbids_tmp);
			}else{
				tid2pdbid.get(tid).add(pdbid);
				tid2pdbid_tmp.get(tid).add(pdbid);
			}
		}
		debug_println("tid2pdbid() "+ (tid2pdbid.size()));
		outputMapAsFile(tid2pdbid, "tid2pdbid");
		/*for (String k: tid2pdbid.keySet()){
			System.out.println(k +" "+ tid2pdbid.get(k));
			
		}*/
		
		PreparedStatement getMod = m_conn.prepareStatement("SELECT modid"
				+ "	FROM ting.isoform_mod"
				+ " WHERE geneid = ? ");
		
		Map<String, Set<String>> m_Geneid2knownsMap = GeneAnnot.getGeneid2knowns(m_conn);
		Set<String> genesym_set = new HashSet<String>();
		List<String> lines2 = FileUtil.getFileAsList("out_interaction_features_3.log");
		BufferedWriter bw1 = new BufferedWriter(new FileWriter("StructureCase.out"));
		bw1.write("No.\tTid1\tGeneId1\tGeneSym1\tProteinId1\t" + 
		          "Tid2\tGeneId2\tGeneSym2\tProteinId2\t" +
		          "Score\n");
		int i = 0;
		int no = 0;
		for (String line2: lines2) {
			i++;
			if (i % 1000 == 0){
				System.out.println(i);
			}
			String[] items = line2.split("\t");
			String tid1 = items[0];
			String geneid1 = items[1];
			String genesym1 = items[2];
			String proteinid1 = items[3];
			String tid2 = items[4];
			String geneid2 = items[5];
			String genesym2 = items[6];
			String proteinid2 = items[7];
			String score = items[8];
			if (!tid2pdbid.keySet().contains(tid1) || !tid2pdbid.keySet().contains(tid2)){
				continue;
			}
			Set<String> known_set = m_Geneid2knownsMap.get(geneid1);
			if (known_set.size() == 1){
				continue;
			}
			
			getMod.setString(1, geneid1);
			ResultSet rs = getMod.executeQuery();
			String mod = "";
			while (rs.next()) {
				mod = rs.getString("modid");
			}
			if (mod.equals("")){
				System.out.println("module not found");
				continue;
			}
			
			Set<String> strucs1_sub = new HashSet<String>();
			strucs1_sub = tid2pdbid.get(tid1);
			Set<String> strucs2_sub = new HashSet<String>();
			strucs2_sub = tid2pdbid.get(tid2);
			Set<String> strucs_sub_tmp = new HashSet<String>();
			strucs_sub_tmp = tid2pdbid_tmp.get(tid1);
			Set<String> strucs_sub_copy = new HashSet<String>();
			for (String c:strucs_sub_tmp){
				strucs_sub_copy.add(c);
			}
			
			strucs_sub_copy.retainAll(strucs2_sub);
			if (strucs_sub_copy.size() > 0){
				strucs1_sub.removeAll(strucs_sub_copy);
				strucs2_sub.removeAll(strucs_sub_copy);
			}
			
			if (strucs1_sub.size() == 0 || strucs2_sub.size() == 0){
				continue;
			}
			
			Set<String> strucs_tmp = new HashSet<String>();
			Set<String> strucs1 = new HashSet<String>();
			Set<String> strucs2 = new HashSet<String>();
			for (String tmp_sub: strucs1_sub){
				String tmp = tmp_sub.substring(0, tmp_sub.indexOf("_"));
				strucs_tmp.add(tmp);
			}
			for (String s1_sub: strucs1_sub){
				String s1 = s1_sub.substring(0, s1_sub.indexOf("_"));
				strucs1.add(s1);
			}
			for (String s2_sub: strucs2_sub){
				String s2 = s2_sub.substring(0, s2_sub.indexOf("_"));
				strucs2.add(s2);
			}

			strucs_tmp.retainAll(strucs2);
			if (strucs_tmp.size() == 0){
				continue;
			}
			
			genesym_set.add(genesym1);
			no++;
			bw1.write(genesym_set.size() +"\t"+ known_set + "\n");
			/*bw1.write(no +"\t"+ tid1 +"\t"+ geneid1 +"\t"+ genesym1 +"\t"+ proteinid1 +"\t"+ struc1_sub_str +"\t"+
					tid2 +"\t"+ geneid2 +"\t"+ genesym2 +"\t"+ proteinid2 +"\t"+ struc2_sub_str +"\t"+
					score + "\n");*/
			for (String s_both: strucs_tmp){

				String struc1_sub_str = "";
				for (String s1_sub: strucs1_sub){
					if (s1_sub.contains(s_both))
						struc1_sub_str += s1_sub + ",";
				}

				String struc2_sub_str = "";
				for (String s2_sub: strucs2_sub){
					if (s2_sub.contains(s_both))
						struc2_sub_str += s2_sub + ",";
				}

				bw1.write(genesym_set.size() +"\t"+ genesym1 +"\t"+ geneid1 +"\t"+ proteinid1 +"\t"+ tid1 +"\t"+ struc1_sub_str +"\t"+
						            genesym2 +"\t"+ geneid2 +"\t"+ proteinid2 +"\t"+ tid2 +"\t"+ struc2_sub_str +"\t"+
						            score + "\n");
				/*bw1.write("\t"+ s_both +"\t"+ tid1 +"\t"+ struc1_sub_str +"\t"
											+ tid2 +"\t"+ struc2_sub_str +"\t" + "\n");*/
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
