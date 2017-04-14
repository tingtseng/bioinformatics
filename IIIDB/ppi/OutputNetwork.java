package edu.nchu.syslab.ppi;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.ResourceBundle;
import java.util.Set;

import edu.nchu.syslab.TextUtil;
import edu.usc.zhoulab.common.MathUtil;
import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.Enrich;
import edu.usc.zhoulab.common.annot.GeneAnnot;
import edu.usc.zhoulab.common.annot.ProteinAnnot;

public class OutputNetwork { 
	static ResourceBundle bundler = ResourceBundle.getBundle("common");
	static String DB_URL = bundler.getString("db.url");
	static String DB_USER = bundler.getString("db.user");
	static String DB_PASS = bundler.getString("db.passwd");
	
	static String hg18_to_uniprot =	"/home/jimliu/ProcessUniprot/2011Mar/hg18_to_uniprot.txt";
	static Connection m_conn = null;
	static boolean m_debug_flag = true;
	static String m_homologene_data = "/home/jimliu/gene/HomoloGene.Oct2010/homologene.data";
	static String m_option;
	
	public OutputNetwork() 
	throws Exception {
		if (m_option.contains("annot")) {
			homology();
			return;
		}		
		if (m_option.contains("eggNOG")) {
			eggNOG();
		}
		if (m_option.contains("studyPfam")) {
			studyPfam(20);
			studyPfam(10);
			studyPfam(5);
			studyPfam(3);
			studyPfam(2);
		}
	}
	
	public static void main(String[] args) 
	throws Exception
	{
		if (args.length != 1) {
			System.out.println("Usage: OutputNetwork option ");
			return;
		}
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);	
		m_option = args[0];
		System.out.println("m_option "+m_option);
		new OutputNetwork();
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


	Map<String, Set<String>> m_goid2tid = null;
	static Set<String> m_intActknowns = null;
	static Set<String> m_geneid_ppi = null;
	static Map<String, Set<String>> m_knowns2partners = null;
	static void addKnownInt (String t1, String t2) {
		Set<String> partners = m_knowns2partners.get(t1);
		if (partners == null) {
			partners = new HashSet<String>();
			m_knowns2partners.put(t1, partners);
		}
		partners.add(t2);		
	}
	
	static public Set<String> getIntActIso(Connection conn) 
	throws Exception{
		m_knowns2partners = new HashMap<String, Set<String>>();
		m_intActknowns = new HashSet<String>();
		m_geneid_ppi = new HashSet<String>();
		Set<String> knownInt = new HashSet<String>();
		Set<String> proteinInt = new HashSet<String>();
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		Map<String, String> known2geneid = GeneAnnot.getKnown2geneid(conn);
		
		List<String> lines = FileUtil.getFileAsList(
				GeneAnnot.getGeneDir() + "intact.Jun2010/intact.txt");
		Set<String> proteins1 = new HashSet<String>();		
		Set<String> proteins2 = new HashSet<String>();		
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

			Set<String> knowns1 = protein2knowns.get(protein1);
			//debug_println("  protein1 "+ protein1 +" -> knowns "+ knowns1);
			if (knowns1 == null && protein1.contains("-")) {
				protein1 = protein1.substring(0, protein1.indexOf('-'));
				knowns1 = protein2knowns.get(protein1);
				//debug_println("    "+ protein1 +" -> knowns "+ knowns1);
			}
			if (knowns1 == null) {
				continue;
			}
			Set<String> knowns2 = protein2knowns.get(protein2);
			//debug_println("  protein2 "+ protein2 +" -> knowns "+ knowns2);
			if (knowns2 == null && protein2.contains("-")) {
				protein2 = protein2.substring(0, protein2.indexOf('-'));
				knowns2 = protein2knowns.get(protein2);
			}
			if (knowns2 == null) {
				continue;
			}			
			proteins1.add(protein1);
			proteins2.add(protein2);
			m_intActknowns.addAll(knowns1);
			m_intActknowns.addAll(knowns2);			
			for(String t1: knowns1) {
				String geneid1 = known2geneid.get(t1);
				for(String t2: knowns2) {
					String geneid2 = known2geneid.get(t2);
					if (geneid1 != null && geneid2 != null) {
						m_geneid_ppi.add(geneid1+":"+geneid2);
						m_geneid_ppi.add(geneid2+":"+geneid1);
					}
					knownInt.add(t1+":"+t2);
					knownInt.add(t2+":"+t1);
					addKnownInt(t1, t2);
					addKnownInt(t2, t1);
				}
			}			
		}
		debug_println("  proteins 1 "+ proteins1.size());
		debug_println("  proteins 2 "+ proteins2.size());
		proteins1.addAll(proteins2);
		debug_println("  proteins merge "+ proteins1.size());		
		debug_println("knownInt "+ (knownInt.size()/2));	
		debug_println("proteinInt "+ (proteinInt.size()/2));	
		debug_println("m_geneid_ppi "+ m_geneid_ppi.size());
		debug_println("m_knowns2partners "+ m_knowns2partners.size());
		checkKnowns2partners(m_knowns2partners);
		return knownInt;
		
	}

	static void checkKnowns2partners(Map<String, Set<String>> m_knowns2partners)
			throws Exception {
		// Map<String, String> known2geneid = GeneAnnot.getKnown2geneid(m_conn);
		Map<String, Set<String>> geneid2knowns = GeneAnnot
				.getGeneid2knowns(m_conn);
		Set<String> sel_gene1 = new HashSet<String>();
		Set<String> sel_gene2 = new HashSet<String>();
		Set<String> sel_trans1 = new HashSet<String>();
		Set<String> sel_trans2 = new HashSet<String>();
		Set<String> sel_trans3 = new HashSet<String>();

		for (String geneid : geneid2knowns.keySet()) {
			Set<String> knowns = new HashSet<String>(geneid2knowns.get(geneid));
			Set<String> annot_tids = new HashSet<String>();
			for (String tid : knowns) {
				Set<String> partners = m_knowns2partners.get(tid);
				if (partners != null) {
					annot_tids.add(tid);
				}
			}
			if (annot_tids.size() < 2) {
				continue;
			}
			sel_gene1.add(geneid);
			boolean found = false;
			Set<String> go1 = null;
			Set<String> go2 = null;
			for (String t : annot_tids)
			// for (String t: knowns)
			{
				if (go1 == null) {
					go1 = m_knowns2partners.get(t);
					if (go1 == null)
						go1 = new HashSet<String>();
				} else {
					if (m_knowns2partners.get(t) == null) {
						go2 = new HashSet<String>();
					} else {
						go2 = new HashSet<String>(m_knowns2partners.get(t));
					}
					if (Math.abs(go1.size() - go2.size()) >= 1) {
						found = true;
					}
					go2.removeAll(go1);
					if (go2.size() >= 1) {
						found = true;
					}
				}
			}
			if (found) {
				sel_gene2.add(geneid);
				sel_trans2.addAll(knowns);
				sel_trans3.addAll(annot_tids);
			}
		}
		FileUtil.outputCollectionAsFile(sel_gene1, "sel_gene1."
				+ sel_gene1.size() + ".log");
		FileUtil.outputCollectionAsFile(sel_gene2, "sel_gene2."
				+ sel_gene2.size() + ".log");
		FileUtil.outputCollectionAsFile(sel_trans1, "sel_trans1."
				+ sel_trans1.size() + ".log");
		FileUtil.outputCollectionAsFile(sel_trans2, "sel_trans2."
				+ sel_trans2.size() + ".log");
		FileUtil.outputCollectionAsFile(sel_trans3, "sel_trans3."
				+ sel_trans3.size() + ".log");
	}

	Map<String, Set<String>> pfam2tids = null;
	Map<String, Set<String>> getTid2Pfam() throws Exception {
		pfam2tids = new HashMap<String, Set<String>>();
		Map<String, Set<String>> tid2Pfam = new HashMap<String, Set<String>>();
		PreparedStatement getTid = m_conn
			.prepareStatement("SELECT distinct name, value"
				+ "	FROM ucsc_hg18.knownToPfam");		
		ResultSet rs = getTid.executeQuery();			
		while (rs.next()) {
			String tid = rs.getString("name");
			String Pfam = rs.getString("value");
			if (tid != null)
			{
				Set<String> Pfam_set = tid2Pfam.get(tid);
				if (Pfam_set == null){
					Pfam_set = new HashSet<String>();
				}
				Pfam_set.add(Pfam);
				tid2Pfam.put(tid, Pfam_set);
			}
			if (Pfam != null)
			{
				Set<String> tid_set = pfam2tids.get(Pfam);
				if (tid_set == null){
					tid_set = new HashSet<String>();
				}
				tid_set.add(tid);
				pfam2tids.put(Pfam, tid_set);
			}
			
		}		
		getTid.close();
		debug_println("getTid2Pfam() tid2Pfam "+ tid2Pfam.size());
		debug_println("  pfam2tids "+ pfam2tids.size());
		return tid2Pfam;
	}
	
	void studyPfam(int freq) throws Exception {
		List<String> tid2pfam_text = new ArrayList<String>();
		debug_println("studyPfam() freq "+freq);
		Map<String, Set<String>> tid2Pfam = getTid2Pfam();
		List<String> pfam1 = new ArrayList<String>();
		for (String p: pfam2tids.keySet()) {
			Set<String> tid_set = pfam2tids.get(p);
			//debug_println("    "+p+"\t"+ tid_set.size());
			if (tid_set.size() > freq) {
				pfam1.add(p);
				for (String t: tid_set) {
					tid2pfam_text.add(t +"\t"+ p);
				}
			}
		}
		debug_println("  pfam1 "+pfam1.size());		
		Set<String> pfam2 = new HashSet<String>();
		Set<String> pfam3 = new HashSet<String>();
		/*
		for (String p1: pfam1) {
			Set<String> tid_set1 = pfam2tids.get(p1);
			for (String p2: pfam1) {
				if (p2.equals(p1)) {
					continue;
				}
				Set<String> tid_set2 = new HashSet<String>(pfam2tids.get(p2));
				tid_set2.retainAll(tid_set1);
				if (tid_set2.size() > freq && 
						!pfam2.contains(p2 +"-"+ p1) && 
						!pfam2.contains(p1 +"-"+ p2)) {					
					String pair = p1 +"-"+p2; 
					pfam2.add(pair);
					for (String t: tid_set2) {
						tid2pfam_text.add(t +"\t"+ pair);
					}
				}
			}			
		}	*/	
		
		for (int i=0; i<pfam1.size(); i++) {
			String p1 = pfam1.get(i);
			Set<String> tid_set1 = pfam2tids.get(p1);
			for (int j=i+1; j<pfam1.size(); j++) {
				String p2 = pfam1.get(j);
				Set<String> tid_set2 = new HashSet<String>(pfam2tids.get(p2));
				tid_set2.retainAll(tid_set1);
				if (tid_set2.size() <= freq) {
					continue;
				}
				pfam2.add(p1 +"-"+p2);	
				for (String t: tid_set2) {
					tid2pfam_text.add(t +"\t"+ p1 +"-"+p2);
				}
				/*
				for (int k=0; k<pfam1.size();k++) {
					String p3 = pfam1.get(k);
					Set<String> tid_set3 = new HashSet<String>(pfam2tids.get(p3));
					tid_set3.retainAll(tid_set2);
					if (tid_set3.size() <= freq) {
						continue;
					}
					pfam3.add(p1 +"-"+p2+"-"+p3);
					for (String t: tid_set3) {
						tid2pfam_text.add(t +"\t"+ p1 +"-"+p2+"-"+p3);
					}				
				}*/							
			}			
		}
		debug_println("  pfam2 "+pfam2.size());		
		debug_println("  pfam3 "+pfam3.size());		
		debug_println("  total domains "+(pfam1.size() + pfam2.size() + pfam3.size()));
		FileUtil.outputCollectionAsFile(tid2pfam_text, "tid2pfam_d2_f"+freq+".txt");
	}
	
	void eggNOG() throws Exception {
		//List<String> annot = FileUtil.getFileAsList(m_option +"/protein.aliases.v2.txt");
		List<String> annot = FileUtil.getFileAsList(m_option +"/accToTaxon.txt");		
		debug_println("eggNOG()  annot "+annot.size());
		Set<String> yeast_proteins = new HashSet<String>();
		Set<String> human_proteins = new HashSet<String>();
		for (String line: annot) {
			if (line.startsWith("#")) {
				continue;
			}
			String[] items = line.split("\t");
			int tax = Integer.parseInt(items[1]);
			String protein = items[0];
			if (tax == GeneAnnot.TAX_YEAST) {
				yeast_proteins.add(protein);
			}
			if (tax == GeneAnnot.TAX_HUMAN) {
				human_proteins.add(protein);
			}			
		}
		annot = null;
		Map<String, Set<String>> yeast2human = 
			getYeast2human(yeast_proteins, human_proteins);
		//yeast2humanIntAct(m_conn, yeast2human);
		debug_println("  yeast_proteins "+yeast_proteins.size());
		debug_println("  human_proteins "+human_proteins.size());
	}

	Map<String, Set<String>> getProtein2nog(Set<String> proteins)
		throws Exception {
		Map<String, Set<String>> protein2nogs = new HashMap<String, Set<String>>();
		List<String> UniProtAC2eggNOG = FileUtil.getFileAsList(m_option +"/UniProtAC2eggNOG.tsv");
		Set<String> nogs = new HashSet<String>();
		for (String line: UniProtAC2eggNOG) {
			String[] items = line.split("\t");
			String p = items[0];
			String[] nog = items[1].split(" ");
			//debug_println("    nog " + nog.length +" "+items[1]);
			for (String n: nog) {
				n = n.replace("fu", "");
				nogs.add(n);
				if (proteins.contains(p)) {
					Set<String> p_nogs = protein2nogs.get(p);
					if (p_nogs == null) {
						p_nogs = new HashSet<String>();
						protein2nogs.put(p, p_nogs);
					}
					p_nogs.add(n);
				}				
			}
		}
		debug_println("  nogs " + nogs.size());		
		debug_println("  protein2nogs " + protein2nogs.size());	
		return protein2nogs;
	}

	Map<String, Set<String>> getNog2proteins(Set<String> proteins)
			throws Exception {
		Map<String, Set<String>> nog2proteins = new HashMap<String, Set<String>>();
		List<String> UniProtAC2eggNOG = FileUtil.getFileAsList(m_option
				+ "/UniProtAC2eggNOG.tsv");
		//Set<String> nogs = new HashSet<String>();
		for (String line : UniProtAC2eggNOG) {
			String[] items = line.split("\t");
			String p = items[0];
			String[] nog = items[1].split(" ");
			// debug_println("    nog " + nog.length +" "+items[1]);
			for (String n : nog) {
				//nogs.add(n);
				n = n.replace("fu", "");
				if (proteins.contains(p)) {
					Set<String> p1 = nog2proteins.get(n);
					if (p1 == null) {
						p1 = new HashSet<String>();
						nog2proteins.put(n, p1);
					}
					p1.add(p);
				}
			}
		}
		//debug_println("  nogs " + nogs.size());
		debug_println("  nog2proteins " + nog2proteins.size());
		return nog2proteins;
	}
	
	Map<String, Set<String>> getYeast2human(Set<String> yeast_proteins, Set<String> human_proteins)
		throws Exception {
		Map<String, Set<String>> yeast2nogs = getProtein2nog(yeast_proteins);
		Map<String, Set<String>> nog2human = getNog2proteins(human_proteins);
		Map<String, Set<String>> yeast2human = new HashMap<String, Set<String>>();
		for (String yeast: yeast2nogs.keySet()) {
			Set<String> human = new HashSet<String>();
			Set<String> nogs = yeast2nogs.get(yeast);
			debug_println("  yeast " + yeast +" -> "+nogs);
			for (String n: nogs) {
				Set<String> p = nog2human.get(n);
				debug_println("    nog " + n +" -> "+p);
				if (p != null) {
					human.addAll(p);
				}
			}
			if (human.size() > 0) {
				yeast2human.put(yeast, human);
			}			
		}
		debug_println("  yeast2nogs " + yeast2nogs.size());		
		debug_println("  nog2human " + nog2human.size());		
		debug_println("  yeast2human " + yeast2human.size());		
		return yeast2human;		
	}
	
	public Set<String> yeast2humanIntAct(Connection conn, Map<String, Set<Hit>> yeast2human) 
	throws Exception{
		m_knowns2partners = new HashMap<String, Set<String>>();
		m_intActknowns = new HashSet<String>();
		m_geneid_ppi = new HashSet<String>();
		Set<String> knownInt = new HashSet<String>();
		Set<String> proteinInt = new HashSet<String>();
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		
		List<String> lines = FileUtil.getFileAsList(
				"/home/ting/IntAct/2010_07_23/psimitab/intact.txt");
				// GeneAnnot.getGeneDir() + "intact.Jun2010/intact.txt");
		Set<String> all_proteins = new HashSet<String>();		
		for (String line: lines) {
			//if (line.startsWith("#") || !line.contains("taxid:4932(Baker's yeast)")) 
			if (line.startsWith("#") || line.contains("taxid:9606(Human)"))
			{
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
			all_proteins.add(protein1);
			all_proteins.add(protein2);
			
			Set<Hit> h_proteins1 = yeast2human.get(protein1);
			Set<Hit> h_proteins2 = yeast2human.get(protein2);
			if (h_proteins1 == null && protein1.contains("-")) {
				protein1 = protein1.substring(0, protein1.indexOf('-'));
				h_proteins1 = yeast2human.get(protein1);
			}
			if (h_proteins2 == null && protein2.contains("-")) {
				protein2 = protein1.substring(0, protein2.indexOf('-'));
				h_proteins2 = yeast2human.get(protein2);
			}
			if (h_proteins1 == null || h_proteins2 == null) {
				continue;
			}
			Set<Hit> knowns1 = new HashSet<Hit>();
			Set<Hit> knowns2 = new HashSet<Hit>();
			for (Hit h1: h_proteins1) {
				Set<String> k1 = protein2knowns.get(h1.protein);
				if (k1 == null) {
					continue;
				}
				for (String known_k1: k1) {
					knowns1.add(new Hit(known_k1, h1.sim));
				}
			}
			for (Hit h2: h_proteins2) {
				Set<String> k2 = protein2knowns.get(h2.protein);
				if (k2 == null) {
					continue;
				}
				for (String known_k2: k2) {
					knowns2.add(new Hit(known_k2, h2.sim));
				}
			}
			
			for(Hit t1: knowns1) {
				for(Hit t2: knowns2) {
					double sim = (t1.sim + t2.sim)/2;
					knownInt.add(t1.protein+"\t"+t2.protein +"\t"+ sim);
					knownInt.add(t2.protein+"\t"+t1.protein +"\t"+ sim);
				}
			}			
		}
		debug_println("Yeast proteins "+ all_proteins.size());		
		debug_println("Yeast proteinInt "+ (proteinInt.size()/2));	
		debug_println("knownInt "+ (knownInt.size()/2));	
		debug_println("m_geneid_ppi "+ m_geneid_ppi.size());
		debug_println("m_knowns2partners "+ m_knowns2partners.size());
		FileUtil.outputCollectionAsFile(knownInt, "knownInt.txt");
		return knownInt;
		
	}

	Map<String, Set<String>> getGeneid2homo(int taxon) throws Exception {
		Map<String, Set<String>> geneid2homo = new HashMap<String, Set<String>>();
		List<String> data = FileUtil.getFileAsList(m_homologene_data);
		Set<String> all_homos = new HashSet<String>();
		for (String line : data) {
			String[] items = line.split("\t");
			int tax = Integer.parseInt(items[1]);
			if (tax != taxon) {
				continue;
			}
			String homo = items[0];
			all_homos.add(homo);
			String geneid = items[2];
			Set<String> homo1 = geneid2homo.get(geneid);
			if (homo1 == null) {
				homo1 = new HashSet<String>();
				geneid2homo.put(geneid, homo1);
			}
			homo1.add(homo);
		}
		debug_println("  all_homos " + all_homos.size());
		debug_println("  geneid2homo " + geneid2homo.size());
		return geneid2homo;
	}	
	
	class Hit {
		String protein;
		double sim;
		public Hit(String protein1, double sim1) {
			protein = protein1;
			sim = sim1;
		}
	}
	void homology() throws Exception {
		double sim_th = 0.3;
		Map<String, Set<Hit>> protein2human = new HashMap<String, Set<Hit>>();
		List<String> mapping = FileUtil.getFileAsList(m_option);		
		debug_println("homology() mapping "+mapping.size());
		debug_println("  sim_th "+sim_th);
		int count = 0;
		for (String line: mapping) {
			String[] items = line.split("\t");
			String p1 = items[2];
			String p2 = items[0];
			double sim = Double.parseDouble(items[25]);
			if (sim <= sim_th) {
				continue;
			}
			count++;
			Set<Hit> hits = protein2human.get(p1);
			if (hits == null) {
				hits = new HashSet<Hit>();
				protein2human.put(p1, hits);
			}
			hits.add(new Hit(p2, sim));			
		}
		//yeast2humanIntAct(m_conn, protein2human);
		protein2humanIntAct(m_conn, protein2human);
		debug_println("mapping count "+count);
	}
	

	public Set<String> protein2humanIntAct(Connection conn, Map<String, Set<Hit>> protein2human) 
	throws Exception{
		String intact = "/home/ting/IntAct/2011_03_09/psimitab/intact.txt";
		debug_println("protein2humanIntAct() intact "+ intact);
		m_knowns2partners = new HashMap<String, Set<String>>();
		m_intActknowns = new HashSet<String>();
		m_geneid_ppi = new HashSet<String>();
		Set<String> knownInt = new HashSet<String>();
		Set<String> proteinInt = new HashSet<String>();
		Map<String, Set<String>> protein2knowns = null;
		if (hg18_to_uniprot == null){
			protein2knowns = GeneAnnot.getProtein2knowns(m_conn);
		}else{
			protein2knowns = GeneAnnot.getProtein2knowns(hg18_to_uniprot);
		}
		
		List<String> lines = FileUtil.getFileAsList(intact);
				// "/home/ting/IntAct/2010_07_23/psimitab/intact.txt");				
		Set<String> all_proteins = new HashSet<String>();		
		for (String line: lines) {
			//if (line.startsWith("#") || !line.contains("taxid:4932(Baker's yeast)")) 
			if (line.startsWith("#") || line.contains("taxid:9606(Human)"))
			{
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
			all_proteins.add(protein1);
			all_proteins.add(protein2);
			
			Set<Hit> knowns1 = protein2human.get(protein1);
			Set<Hit> knowns2 = protein2human.get(protein2);
			if (knowns1 == null || knowns2 == null) {
				continue;
			}
			for(Hit t1: knowns1) {
				for(Hit t2: knowns2) {
					double sim = (t1.sim + t2.sim)/2;
					knownInt.add(t1.protein+"\t"+t2.protein +"\t"+ sim);
					knownInt.add(t2.protein+"\t"+t1.protein +"\t"+ sim);
				}
			}			
		}
		debug_println("Yeast proteins "+ all_proteins.size());		
		debug_println("Yeast proteinInt "+ (proteinInt.size()/2));	
		debug_println("knownInt "+ (knownInt.size()/2));	
		debug_println("m_geneid_ppi "+ m_geneid_ppi.size());
		debug_println("m_knowns2partners "+ m_knowns2partners.size());
		FileUtil.outputCollectionAsFile(knownInt, "knownInt.txt");
		return knownInt;
		
	}	
}
