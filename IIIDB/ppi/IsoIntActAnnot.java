package edu.nchu.syslab.ppi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.ResourceBundle;
import java.util.Set;

import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.Enrich;
import edu.usc.zhoulab.common.annot.GeneAnnot;

public class IsoIntActAnnot {
	static ResourceBundle bundler = ResourceBundle.getBundle("common");
	static String DB_URL = bundler.getString("db.url");
	static String DB_USER = bundler.getString("db.user");
	static String DB_PASS = bundler.getString("db.passwd");

	static boolean m_debug_flag = true;
	public static boolean USE_SYS_EXIT = true;
	public static String fileName = null;
	static Connection m_conn = null;
	
	Random random = new Random();

	public static void main(String[] args) throws Exception {
		if (args.length != 1) {
			System.out.println("Usage: IsoIntActAnnot fileName");
									// IsoIntActAnnot D:\190\IntAct_webserver\db\isoform_mod.txt
			return;
		}
		String fileName = args[0];
		System.out.println("fileName " + fileName);
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);
		new IsoIntActAnnot();
	}
	
	public IsoIntActAnnot() throws Exception {
		//OutModule2Tid();
		//EnrichPath();
		//EnrichGO();
		//ProcessModuleFile2();
		//ProcessCocitation();
		//OutModuleList();
		//EnrichComplex();
		
		//OutModule2Tid_new();
		//EnrichPath_new();
		//EnrichGO_new();
		//EnrichComplex_new();
		
		//ProcessModuleFile2_new();
		
		//OutRandomModule_new();
		
		//OutModuleList_new();
		
		//getFileSize();
		
		OutModuleList();
		
	}
	
	void ProcessCocitation() throws Exception {
		List<String> lines = FileUtil.getFileAsList("gene2pubmed");
		//9	1246500	9873079
		Map<String, Set<String>> geneid2pubmedid = new HashMap<String, Set<String>>();
		int j = 0;
		for (String line: lines) {
			j++;
			if (j % 1000000 == 0)
				System.out.println(j);
			String[] items = line.split("\t");
			String geneid = items[1];
			String pubmedid = items[2];
			Set<String> pubmedids = geneid2pubmedid.get(geneid);
			if (pubmedids == null) {
				pubmedids = new HashSet<String>();
				geneid2pubmedid.put(geneid, pubmedids);
			}
			pubmedids.add(pubmedid);
		}
		debug_println("geneid2pubmedid.size() "+ (geneid2pubmedid.size()));
		List<String> lines2 = FileUtil.getFileAsList("out_interaction_features.log");
		BufferedWriter bw = new BufferedWriter(new FileWriter("co_citation"));
		int i = 0;
		for (String line: lines2) {
			i++;
			String[] items = line.split("\t");
			String geneid1 = items[1];
			String geneid2 = items[5];
			Set<String> citations1 = geneid2pubmedid.get(geneid1);
			Set<String> citations2 = geneid2pubmedid.get(geneid2);
			if (citations1 == null || citations2 == null){
				System.out.println("null");
				if (i % 1000 == 0)
					System.out.println(i);
				continue;
			}
			if (citations1.retainAll(citations2)){
				bw.write(geneid1 + "\t" + geneid2 + "\t" + citations1 + "\n");
				System.out.print(geneid1 + "\t" + geneid2 + "\t" + citations1 + "\n");
			}else{
				if (i % 1000 == 0)
					System.out.println(i);
			}
		}
		bw.close();
	}
	
	void EnrichPath() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_kegg("hsa/hsa_pathway.list");
		BufferedWriter bw_path = new BufferedWriter(new FileWriter("out_path_annot.log"));
		List<String> lines = FileUtil.getFileAsList("isoform_mod.txt");
		for (String line: lines) {
			StringBuffer output = new StringBuffer();
			String[] items = line.split("\t");
			String mod = items[0];
			int gene_num = Integer.parseInt(items[2]);
			if (gene_num > 12){
				continue;
			}
			String tid_string = items[4];
			String[] tids = tid_string.split(", ");
			Set<String> geneids = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				geneids.add(geneid);
			}
			Set<String> path_Result = Enrich.enrich_path(geneids, output, 0.001, 500, 2);
			//hsResult.add(ge.m_id + " " + name+ "\t" + df1.format(ge.getP()) + "\t" + GeneAnnot.getGeneSym(ge.m_local_geneid));
			if (path_Result.size() == 0)
				continue;
			for (String result: path_Result){
				String[] results = result.split("\t");
				String path_id_name = results[0];
				String pvalue = results[1];
				String associ_genes = results[2];
				associ_genes = associ_genes.replace("[", "");
				associ_genes = associ_genes.replace("]", "");
				bw_path.write("mod-" + mod + "\t" + path_id_name + "\t" + associ_genes + "\t" + pvalue + "\n");
				//debug_println("mod-" + mod + "\t" + path_id_name + "\t" + pvalue + "\t" + associ_genes);
			}	
		}	
		bw_path.close();
	}

	void EnrichGO() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true, "go2gene.9606");
		BufferedWriter bw_go = new BufferedWriter(new FileWriter("out_GO_annot.log"));
		List<String> lines = FileUtil.getFileAsList("isoform_mod.txt");
		for (String line: lines) {
			StringBuffer output = new StringBuffer();
			String[] items = line.split("\t");
			String mod = items[0];
			int gene_num = Integer.parseInt(items[2]);
			if (gene_num > 12){
				continue;
			}
			String tid_string = items[4];
			String[] tids = tid_string.split(", ");
			Set<String> geneids = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				geneids.add(geneid);
			}
			
			Set<String> GO_Result = Enrich.enrich(geneids, output, 0.001, 500, 2, "GO");
			if (GO_Result.size() == 0)
				continue;
			for (String result: GO_Result){
				//result.add(ge.m_id +"\t"+df1.format(ge.getP())+"\t" + GeneAnnot.getGeneSym(ge.m_local_geneid));
				String[] results = result.split("\t");
				String go_id = results[0];
				String pvalue = results[1];
				String associ_genes = results[2];
				associ_genes = associ_genes.replace("[", "");
				associ_genes = associ_genes.replace("]", "");
				bw_go.write("mod-" + mod + "\t" + go_id + "\t" + associ_genes + "\t" + pvalue + "\n");
				//debug_println("mod-" + mod + "\t" + go_id + "\t" + pvalue + "\t" + associ_genes);
			}
		}
		bw_go.close();
	}

	void EnrichPath_new() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_kegg("hsa/hsa_pathway.list");
		BufferedWriter bw_path = new BufferedWriter(new FileWriter("out_path_annot.log"));
		List<String> lines = FileUtil.getFileAsList("out_d0.7_random.txt");
		int mod = 0;
		for (String line: lines) {
			mod++;
			StringBuffer output = new StringBuffer();
			String[] tids = line.split(" ");
			Set<String> geneids = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				geneids.add(geneid);
			}
			Set<String> path_Result = Enrich.enrich_path(geneids, output, 0.01, 500, 2);
			//hsResult.add(ge.m_id + " " + name+ "\t" + df1.format(ge.getP()) + "\t" + GeneAnnot.getGeneSym(ge.m_local_geneid));
			if (path_Result.size() == 0)
				continue;
			for (String result: path_Result){
				System.out.println(result);
				String[] results = result.split("\t");
				String path_id_name = results[0];
				String pvalue = results[1];
				String associ_genes = results[2];
				associ_genes = associ_genes.replace("[", "");
				associ_genes = associ_genes.replace("]", "");
				bw_path.write("mod-" + mod + "\t" + path_id_name + "\t" + associ_genes + "\t" + pvalue + "\n");
				//debug_println("mod-" + mod + "\t" + path_id_name + "\t" + pvalue + "\t" + associ_genes);
			}	
		}	
		bw_path.close();
	}
	
	void EnrichGO_new() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true, "go2gene.9606");
		BufferedWriter bw_go = new BufferedWriter(new FileWriter("out_GO_annot.log"));
		List<String> lines = FileUtil.getFileAsList("out_d0.7_random.txt");
		int mod = 0;
		for (String line: lines) {
			mod++;
			StringBuffer output = new StringBuffer();
			String[] tids = line.split(" ");
			Set<String> geneids = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				geneids.add(geneid);
			}
			Set<String> GO_Result = Enrich.enrich(geneids, output, 0.001, 500, 2, "GO");
			//result.add(ge.m_id +"\t"+df1.format(ge.getP())+"\t" + GeneAnnot.getGeneSym(ge.m_local_geneid));
			if (GO_Result.size() == 0)
				continue;
			for (String result: GO_Result){	
				System.out.println(result);
				String[] results = result.split("\t");
				String go_id = results[0];
				String pvalue = results[1];
				String associ_genes = results[2];
				associ_genes = associ_genes.replace("[", "");
				associ_genes = associ_genes.replace("]", "");
				bw_go.write("mod-" + mod + "\t" + go_id + "\t" + associ_genes + "\t" + pvalue + "\n");
				//debug_println("mod-" + mod + "\t" + go_id + "\t" + pvalue + "\t" + associ_genes);
			}
		}
		bw_go.close();
	}
	
	void EnrichComplex_new() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true, "go2gene.9606");
		//Map<String, Set<String>> gene2ProtComplex = Enrich.parse_gene2ProtComplex("CORUM", "9606");
		Enrich.parse_gene2ProtComplex("CORUM", "9606");
		//String complexid = source + ":" + items[0].trim() +" "+ items[1];
		BufferedWriter bw_go = new BufferedWriter(new FileWriter("out_Complex_enrich.log"));
		List<String> lines = FileUtil.getFileAsList("out_d0.7_random.txt");
		int mod = 0;
		for (String line: lines) {
			mod++;
			StringBuffer output = new StringBuffer();
			String[] tids = line.split(" ");
			Set<String> geneids = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				geneids.add(geneid);
			}
			Set<String> GO_Result = Enrich.enrich(geneids, output, 0.01, 500, 2, "protein_complex");
			if (GO_Result.size() == 0)
				continue;
			for (String result: GO_Result){
				//result.add(ge.m_id +"\t"+df1.format(ge.getP())+"\t" + GeneAnnot.getGeneSym(ge.m_local_geneid));
				System.out.println(result);
				String[] results = result.split("\t");
				String go_id = results[0];
				String pvalue = results[1];
				String associ_genes = results[2];
				associ_genes = associ_genes.replace("[", "");
				associ_genes = associ_genes.replace("]", "");
				bw_go.write("mod-" + mod + "\t" + go_id + "\t" + associ_genes + "\t" + pvalue + "\n");
			}
		}
		bw_go.close();
	}
	
	void EnrichComplex() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true, "go2gene.9606");
		//Map<String, Set<String>> gene2ProtComplex = Enrich.parse_gene2ProtComplex("CORUM", "9606");
		Enrich.parse_gene2ProtComplex("CORUM", "9606");
		BufferedWriter bw_go = new BufferedWriter(new FileWriter("out_Complex_enrich.log"));
		List<String> lines = FileUtil.getFileAsList("isoform_mod.txt");
		for (String line: lines) {
			StringBuffer output = new StringBuffer();
			String[] items = line.split("\t");
			String mod = items[0];
			int gene_num = Integer.parseInt(items[2]);
			if (gene_num > 12){
				continue;
			}
			String tid_string = items[4];
			String[] tids = tid_string.split(", ");
			Set<String> geneids = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				geneids.add(geneid);
			}
			Set<String> GO_Result = Enrich.enrich(geneids, output, 0.001, 500, 2, "protein_complex");
			if (GO_Result.size() == 0)
				continue;
			for (String result: GO_Result){
				//result.add(ge.m_id +"\t"+df1.format(ge.getP())+"\t" + GeneAnnot.getGeneSym(ge.m_local_geneid));
				/*String[] results = result.split("\t");
				String go_id = results[0];
				String pvalue = results[1];
				String associ_genes = results[2];
				associ_genes = associ_genes.replace("[", "");
				associ_genes = associ_genes.replace("]", "");
				bw_go.write("mod-" + mod + "\t" + go_id + "\t" + associ_genes + "\t" + pvalue + "\n");*/
				//debug_println("mod-" + mod + "\t" + go_id + "\t" + pvalue + "\t" + associ_genes);
				debug_println("mod-" + mod + "\t" + result);
			}
		}
		bw_go.close();
	}
	
	
	void OutModule2Tid() throws Exception {
		//Set<String> mod2genesym = new HashSet<String>();
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		List<String> lines = FileUtil.getFileAsList("isoform_mod.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter("out_mod2tid.log"));
		for (String line: lines) {
			String[] items = line.split("\t");
			String mod = items[0];
			int gene_num = Integer.parseInt(items[2]);
			if (gene_num > 12){
				continue;
			}
			String tid_string = items[4];
			String[] tids = tid_string.split(", ");
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				String tid = items1[1];
				tid = tid.substring(0, tid.length()-1);
				bw.write("mod-" + mod + "\t" + geneid + "\t" + genesym + "\t" + tid + "\n");
				//mod2genesym.add("mod-" + mod + "\t" + genesym);
			}
		}
		bw.close();
		//return mod2genesym;
	}
	
	void OutModule2Tid_new() throws Exception {
		//Set<String> mod2genesym = new HashSet<String>();
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		List<String> lines = FileUtil.getFileAsList("out_d0.7.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter("out_mod2tid.log"));
		int mod = 0;
		for (String line: lines) {
			mod++;
			String[] tids = line.split(" ");
			for (int i = 0; i < tids.length; i++){
				String genesym_tid = tids[i];
				String[] items1 = genesym_tid.split("\\(");
				String genesym = items1[0];
				String geneid = GeneAnnot.getGeneid(genesym);
				String tid = items1[1];
				tid = tid.substring(0, tid.length()-1);
				bw.write("mod-" + mod + "\t" + geneid + "\t" + genesym + "\t" + tid + "\n");
				//mod2genesym.add("mod-" + mod + "\t" + genesym);
			}
		}
		bw.close();
		//return mod2genesym;
	}
	
	void getFileSize() throws Exception {
               
        File f = new File("./graphic//", "mod-1_S.gif");
        if(f.exists()){
            System.out.println("file name¡G" + f.getAbsolutePath());
            System.out.println("file size¡G" + f.length()/1024 + "KB");
        } else {
            f.getParentFile().mkdirs();
            try{
                f.createNewFile();
            } catch(IOException e){
                e.printStackTrace();
            }
        }
    
	}
	
	void ProcessModuleFile2_new() throws Exception {
		List<String> interact_list = FileUtil.getFileAsList("interaction_30.sif");
		// RPS4Y1(uc004fqi.1) pp HNRNPK(uc004ang.2)
		Set<String> interact_set = new HashSet<String>();
		for (String tid2tid : interact_list) {
			String[] items = tid2tid.split(" pp ");
			String tid1 = items[0];
			tid1 = tid1.trim();
			String tid2 = items[1];
			if (!interact_set.contains(tid2 + ":" + tid1)){
				interact_set.add(tid1 + ":" + tid2);
			}
		}
		List<String> lines = FileUtil.getFileAsList("out_d0.7.txt");
		int mod = 0;
		for (String line: lines) {
			mod++;
			System.out.println("mod-" + mod);
			String[] tids = line.split(" ");
			BufferedWriter bw_S = new BufferedWriter(new FileWriter("./graphic//" + mod + "_S.net"));
			bw_S.write("digraph world { \n");
			bw_S.write("size=\"6,6\"; \n");
			bw_S.write("ratio = 1; \n");
			bw_S.write("fixedsize = true; \n");
			bw_S.write("overlap = false; \n");
			bw_S.write("splines = true; \n");
			bw_S.write("edge [dir=none color=blue style=bold]; \n");
			BufferedWriter bw_B = new BufferedWriter(new FileWriter("./graphic//" + mod + "_B.net"));
			bw_B.write("digraph world { \n");
			bw_B.write("size=\"40,40\"; \n");
			bw_B.write("ratio = 1; \n");
			bw_B.write("fixedsize = false; \n");
			bw_B.write("overlap = false; \n");
			bw_B.write("splines = true; \n");
			bw_B.write("edge [dir=none color=blue style=bold]; \n");
			Set<String> mod_gene_set = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				mod_gene_set.add(tids[i]);
				bw_B.write("\""+tids[i]+"\" \n");
				bw_S.write("\""+tids[i]+"\" \n");
			}
			for (String i:interact_set){
				String[] item1s = i.split(":");
				String tid1 = item1s[0];
				String tid2 = item1s[1];
				if (mod_gene_set.contains(tid1) && mod_gene_set.contains(tid2)){
					bw_B.write("\""+tid1+"\" -> \""+tid2+"\" \n");
					bw_S.write("\""+tid1+"\" -> \""+tid2+"\" \n");
				}
			}
			bw_S.write("}");
			bw_B.write("}");
			bw_B.close();
			bw_S.close();
			
			String cmdline_B = "C://Program Files (x86)//Graphviz 2.28//bin//dot.exe " +
					"-Gmaxiter=99999999 -Tpdf -o D://190//IntAct_8//Module//graphic//mod-"+mod+"_B.pdf " +
					"D://190//IntAct_8//Module//graphic//"+mod+"_B.net";
			Process ps_B = Runtime.getRuntime().exec(cmdline_B);
			try {
				ps_B.waitFor();
			} catch (InterruptedException x) {
			}
			
			String cmdline_S = "C://Program Files (x86)//Graphviz 2.28//bin//dot.exe " +
					"-Gmaxiter=99999999 -Tgif -o D://190//IntAct_8//Module//graphic//mod-"+mod+"_S.gif " +
					"D://190//IntAct_8//Module//graphic//"+mod+"_S.net";
			Process ps_S = Runtime.getRuntime().exec(cmdline_S);
			try {
				ps_S.waitFor();
			} catch (InterruptedException x) {
			}
		}
	}
	
	void ProcessModuleFile2() throws Exception {
		List<String> interact_list = FileUtil.getFileAsList("interaction.sif");
		// RPS4Y1(uc004fqi.1) pp HNRNPK(uc004ang.2)
		Set<String> interact_set = new HashSet<String>();
		for (String tid2tid : interact_list) {
			String[] items = tid2tid.split(" pp  ");
			String tid1 = items[0];
			tid1 = tid1.trim();
			String tid2 = items[1];
			if (!interact_set.contains(tid2 + ":" + tid1)){
				interact_set.add(tid1 + ":" + tid2);
			}
		}
		List<String> lines = FileUtil.getFileAsList("isoform_mod.txt");
		for (String line: lines) {
			String[] items = line.split("\t");
			String mod = items[0];
			int gene_num = Integer.parseInt(items[2]);
			if (gene_num > 12){
				continue;
			}
			String tid_string = items[4];
			String[] tids = tid_string.split(", ");
			BufferedWriter bw_S = new BufferedWriter(new FileWriter("D://Module//" + mod + "_S.net"));
			bw_S.write("digraph world { \n");
			bw_S.write("size=\"6,6\"; \n");
			bw_S.write("ratio = 1; \n");
			bw_S.write("fixedsize = true; \n");
			bw_S.write("overlap = false; \n");
			bw_S.write("splines = true; \n");
			bw_S.write("edge [dir=none color=blue style=bold]; \n");
			BufferedWriter bw_B = new BufferedWriter(new FileWriter("D://Module//" + mod + "_B.net"));
			bw_B.write("digraph world { \n");
			bw_B.write("size=\"15,15\"; \n");
			bw_B.write("ratio = 1; \n");
			bw_B.write("fixedsize = true; \n");
			bw_B.write("overlap = false; \n");
			bw_B.write("splines = true; \n");
			bw_B.write("edge [dir=none color=blue style=bold]; \n");
			Set<String> mod_gene_set = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				mod_gene_set.add(tids[i]);
				bw_B.write("\""+tids[i]+"\" \n");
				bw_S.write("\""+tids[i]+"\" \n");
			}
			for (String i:interact_set){
				String[] item1s = i.split(":");
				String tid1 = item1s[0];
				String tid2 = item1s[1];
				if (mod_gene_set.contains(tid1) && mod_gene_set.contains(tid2)){
					bw_B.write("\""+tid1+"\" -> \""+tid2+"\" \n");
					bw_S.write("\""+tid1+"\" -> \""+tid2+"\" \n");
				}
			}
			bw_S.write("}");
			bw_B.write("}");
			bw_B.close();
			bw_S.close();
			
			String cmdline_B = "C://Program Files (x86)//Graphviz 2.28//bin//dot.exe " +
					"-Gmaxiter=99999999 -Tgif -o D://Module//mod-"+mod+"_B.gif D://Module//"+mod+"_B.net";
			Process ps_B = Runtime.getRuntime().exec(cmdline_B);
			try {
				ps_B.waitFor();
			} catch (InterruptedException x) {
			}
			
			String cmdline_S = "C://Program Files (x86)//Graphviz 2.28//bin//dot.exe " +
					"-Gmaxiter=99999999 -Tgif -o D://Module//mod-"+mod+"_S.gif D://Module//"+mod+"_S.net";
			Process ps_S = Runtime.getRuntime().exec(cmdline_S);
			try {
				ps_S.waitFor();
			} catch (InterruptedException x) {
			}
		}
		
		//digraph world { 
		//size="15,15"; 
		//"ALDH2";
		//"ALDH9A1";
		//"ALDH9A1" -> "ALDH3A2" [dir=none color="black"];
		//"CYP26A1";
		//"CYP26A1" -> "ALDH1A2" [dir=none color="black"];
		//"ALDH9A1";
		//} 
	}
	
	void OutModuleList() throws Exception {
		List<String> interact_list = FileUtil.getFileAsList("interaction_30.sif");
		// RPS4Y1(uc004fqi.1) pp HNRNPK(uc004ang.2)
		Set<String> interact_set = new HashSet<String>();
		for (String tid2tid : interact_list) {
			String[] items = tid2tid.split(" pp ");
			String tid1 = items[0];
			tid1 = tid1.trim();
			String tid2 = items[1];
			if (!interact_set.contains(tid2 + ":" + tid1)){
				interact_set.add(tid1 + ":" + tid2);
			}
		}
		List<String> lines = FileUtil.getFileAsList("out_d0.7.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter("ModuleList_30"));
		bw.write("Module\t#Gene\tGene Symbol\t#Transcript\tGene Symbol(Transcript ID)\t#Interaction\tInteractions\n");
		int mod = 0;
		for (String line: lines) {
			mod++;
			bw.write("mod-" +mod+ "\t");
			String[] tids = line.split(" ");
			Set<String> mod_sym_set = new HashSet<String>();
			Set<String> mod_gene_set = new HashSet<String>();
			Set<String> match_intact_set = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String tid = tids[i];
				String[] syms = tid.split("\\(");
				mod_sym_set.add(syms[0]);
				mod_gene_set.add(tid);
			}
			bw.write(mod_sym_set.size() + "\t");
			for (String s: mod_sym_set){
				bw.write(s + " ");
			}
			bw.write("\t" + mod_gene_set.size() + "\t");
			for (String g: mod_gene_set){
				bw.write(g + " ");
			}
			bw.write("\t");
			for (String i:interact_set){
				String[] item1s = i.split(":");
				String tid1 = item1s[0];
				String tid2 = item1s[1];
				if (mod_gene_set.contains(tid1) && mod_gene_set.contains(tid2)){
					match_intact_set.add(tid1+":"+tid2);		
				}
			}
			bw.write(match_intact_set.size() + "\t");
			for (String m: match_intact_set){
				bw.write(m + " ");
			}
			bw.write("\n");
		}
		bw.close();
	}
	
	void OutModuleList_new() throws Exception {
		List<String> lines = FileUtil.getFileAsList("out_d0.7.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter("ModuleList"));
		bw.write("Module\tGene Symbol\tGene Symbol(Transcript ID)\n");
		int mod = 0;
		for (String line: lines) {
			mod++;
			bw.write("mod-" +mod+ "\t");
			String[] tids = line.split(" ");
			Set<String> mod_sym_set = new HashSet<String>();
			Set<String> mod_gene_set = new HashSet<String>();
			for (int i = 0; i < tids.length; i++){
				String tid = tids[i];
				String[] syms = tid.split("\\(");
				mod_sym_set.add(syms[0]);
				mod_gene_set.add(tid);
			}
			for (String s: mod_sym_set){
				bw.write(s + " ");
			}
			bw.write("\t");
			for (String g: mod_gene_set){
				bw.write(g + " ");
			}
			bw.write("\n");
		}
		bw.close();
	}
	
	void OutRandomModule_new() throws Exception {
		GeneAnnot.parse_gene_info("" + GeneAnnot.TAX_HUMAN);
		Enrich.parse_gene2go("" + GeneAnnot.TAX_HUMAN, true);
		Map<String, String> m_Known2geneidsMap = GeneAnnot.getKnown2geneid(m_conn);
		List<String> lines = FileUtil.getFileAsList("out_d0.7.txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter("out_d0.7_random.txt"));
		List<String> tid_list = FileUtil.getFileAsList("sel_tids.log");
		int mod = 0;
		for (String line: lines) {
			//KIAA0776(uc003por.1) UBA5(uc003epa.2) UFM1(uc001uwu.1)
			mod++;
			String[] tids = line.split(" ");
			for (int i = 0; i < tids.length; i++){
				int random_i = random.nextInt(tid_list.size());
				String tid = tid_list.get(random_i);
				String geneid = m_Known2geneidsMap.get(tid);
				String genesym = GeneAnnot.getGeneSym(geneid);
				bw.write(genesym + "(" + tid + ") ");
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
