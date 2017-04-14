package edu.nchu.syslab.ppi;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.Set;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.HypergeometricDistribution;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;
import edu.usc.zhoulab.common.File.FileUtil;
import edu.usc.zhoulab.common.annot.Enrich;
import edu.usc.zhoulab.common.annot.GeneAnnot;

public class EnrichComplex {
	static ResourceBundle bundler = ResourceBundle.getBundle("common");
	static String DB_URL = bundler.getString("db.url");
	static String DB_USER = bundler.getString("db.user");
	static String DB_PASS = bundler.getString("db.passwd");
	static Connection m_conn = null;

	public static String option = null;
	public static String m_mod_file = null;
	static boolean m_debug_flag = true;
	static DecimalFormat m_df1 = new DecimalFormat();
	static Map<String, Set<String>> m_protein2goMap = null; // gene id -> go name HashSet
	static Map<String, Set<String>> m_go2proteinMap = null; // go name -> gene id HashSet
	static Map<String, Double> m_PvalueCache = new HashMap<String, Double>(); 
	// use cache to improve the performance

	public static void init() {
		m_PvalueCache = new HashMap<String, Double>();
	}

	public EnrichComplex() {
		m_df1.setMaximumFractionDigits(2);
	}

	public static void main(String[] args) {
		if (args.length != 2) {
			System.out.println("Usage: Enrich option mod");
			return;
		}
		String option = args[0];
		m_mod_file = args[1];
		// Enrich Enrich1 = new Enrich();
		try {
			Class.forName("com.mysql.jdbc.Driver").newInstance();
			m_conn = DriverManager.getConnection(DB_URL, DB_USER, DB_PASS);		
			
			if (option.equals("parse_complex")) {
				init();
				parse_complex();
				if (m_mod_file.contains("PPI")) {
					enrich_mod_ppi(0.01);
					enrich_mod_ppi(0.001);
					enrich_mod_ppi(0.0001);
					enrich_mod_ppi(0.00001);
					enrich_mod_ppi(0.000001);
				}
				else {
					enrich_mod_iii(0.01);
					enrich_mod_iii(0.001);
					enrich_mod_iii(0.0001);
					enrich_mod_iii(0.00001);
					enrich_mod_iii(0.000001);
					
				}
				// parse_gene2TF();
				// parse_gene2TF_pheno();
				// parse_gene2Encode5K();
				return;
			}
		} catch (Exception t) {
			t.printStackTrace();
			System.exit(1);
		}

	}

	static protected void debug_println(String str) {
		if (m_debug_flag) {
			System.out.println(str);
		}
	}

	static protected void debug_print(String str) {
		if (m_debug_flag) {
			System.out.print(str);
		}
	}

	// enrich analysis
	// genes: a set of input gene IDs
	// output: the text format for enrichment result
	// p_th: p-value threshold for enrichment analysis
	// max_go_size: skip the GOs with gene nubmer >= max_go_size
	static public Set<String> enrich(Set<String> genes1, StringBuffer output,
			double p_th, int max_go_size, int min_local_size) {
		return enrich(genes1, output, p_th, max_go_size, min_local_size, null,
				1, "GO");
	}

	static public Set<String> enrich(Set<String> genes1, StringBuffer output,
			double p_th, int max_go_size, int min_local_size, String type) {
		return enrich(genes1, output, p_th, max_go_size, min_local_size, null,
				1, type);
	}

	// num_phenotypes is used for FDR adjustment
	// type is TF, GO, or miRNA
	static public Set<String> enrich(Set<String> genes1, StringBuffer output,
			double p_th, int max_go_size, int min_local_size,
			Map<String, Double> go2pvalue, int num_phenotypes, String type) {
		Set<String> result = new HashSet<String>();
		Set<String> genes = new HashSet<String>(genes1);
		genes.retainAll(m_protein2goMap.keySet());
		DecimalFormat df1 = new DecimalFormat();
		df1.applyPattern("0.0000E0");
		Map<String, GoEntry> hmGoEntry = new HashMap<String, GoEntry>();
		for (String go_id : m_go2proteinMap.keySet()) {
			Set<String> global = m_go2proteinMap.get(go_id); // all of genes
															// belong to this
															// go_id term.
			int global_count = global.size();
			if (global_count >= max_go_size) {
				continue;
			}
			Set<String> local = new HashSet<String>(global);
			local.retainAll(genes); // interection genes
			if (local.size() < min_local_size) {
				continue;
			}
			GoEntry ge = new GoEntry();
			ge.m_id = go_id;
			ge.m_global_genes = global_count;
			ge.setLocalprotein(local);
			hmGoEntry.put(go_id, ge);
		}

		if (hmGoEntry.size() == 0) {
			return result;
		}
		/*
		 * output.append("#original input genes " + genes.size() +
		 * ", #total genes " + m_protein2goMap.size() +
		 * ", #max genes of phenotype " + max_go_size + ", #min intersection " +
		 * min_local_size + ", #phenotypes " + hmGoEntry.size() + "\n");
		 */
		// debug_println("hm="+hm);
		Collection<GoEntry> col = hmGoEntry.values();
		GoEntry[] goArray = new GoEntry[col.size()];
		col.toArray(goArray);
		for (int i = 0; i < goArray.length; i++) {
			GoEntry ge = goArray[i];
			String key = ge.m_global_genes + ";" + genes.size() + ";"
					+ ge.m_local_protein.size() + ";";
			if (m_PvalueCache.containsKey(key)) {
				ge.m_p_value = ((Double) m_PvalueCache.get(key)).doubleValue();
			} else {
				HypergeometricDistribution hyper = new HypergeometricDistributionImpl(
						m_protein2goMap.size(), ge.m_global_genes, genes.size());
				try {
					ge.m_p_value = 1.0 - hyper
							.cumulativeProbability(ge.m_local_protein.size() - 1);
					if (ge.m_p_value < 0) {
						ge.m_p_value = 0;
					}
					m_PvalueCache.put(key, new Double(ge.m_p_value));
				} catch (MathException ex) {
					ex.printStackTrace();
					System.exit(1);
				}
			}
			ge.m_p_value *= num_phenotypes;
		}
		sortGoEntry(goArray);
		output.append("Rank" + "\t" +  "Factor Name" + "\t" + "Total" + "\t" + "Found" + "\t" + "P-value" + "\t"
				+ "Transcription Factor Binding Sites (Gene IDs)" + "\t"+ "Transcription Factor Binding Sites (Gene Symbols)" + "\n");

		for (int i = 0; i < goArray.length; i++) {
			GoEntry ge = goArray[i];
			if (go2pvalue != null) {
				go2pvalue.put(ge.m_id, ge.getP());
			}
			
			String gene_list = "";
			
			if (ge.getP() < p_th && ge.m_local_protein.size() >= min_local_size) {				
				result.add(ge.m_id + " (" + df1.format(ge.getP()) + ")");
				output.append((i + 1) + "\t"  + ge.m_id
						+ "\t" + ge.m_global_genes + "\t"
						+ ge.m_local_protein.size() + "\t"
						+ df1.format(ge.getP()) + "\t" 
						+ ge.m_local_protein.toString().replace('[',' ').replace(']',' ')+ "\t" 
						+ gene_list.replace('[',' ').replace(']',' ')
						+ "\n");
			}
		}		
		return result;
	}

	static class GoEntry {
		String m_id = null; // key
		int m_global_genes = 0;
		double m_p_value = 0.0;
		public Set<String> m_local_protein = new HashSet<String>(); // for gene
																	// ID

		public double getP() {
			return m_p_value;
		}

		public void addLocalprotein(String id) {
			m_local_protein.add(id);
		}

		public void setLocalprotein(Set<String> ids) {
			m_local_protein = ids;
		}
	}

	// sort GO entries
	public static void sortGoEntry(GoEntry[] data) {
		Arrays.sort(data, new Comparator<GoEntry>() {
			public int compare(GoEntry o1, GoEntry o2) {
				double diff = ((GoEntry) o1).m_p_value
						- ((GoEntry) o2).m_p_value;
				if (diff > 0) {
					return 1;
				} else if (diff < 0) {
					return -1;
				} else {
					return 0;
				}
			}
		});
	}

	static public Map<String, Set<String>> parse_complex() throws IOException {
		init();
		m_protein2goMap = new HashMap<String, Set<String>>(); // protein -> GO HashSet
		m_go2proteinMap = new HashMap<String, Set<String>>(); // go_id -> protein HashSet
		BufferedReader inputFile = new BufferedReader(new FileReader(
				"/home/ting/IntAct/Complex/pid2complex"));
		String line;
		while ((line = inputFile.readLine()) != null) {
			String[] items = line.split("\t");
			String protein = items[0];
			String complex = items[1].trim();

			Set<String> complexs = m_protein2goMap.get(protein);
			if (complexs == null) {
				complexs = new HashSet<String>();
				m_protein2goMap.put(protein, complexs);
			}
			complexs.add(complex);

			Set<String> proteins = m_go2proteinMap.get(complex);
			if (proteins == null) {
				proteins = new HashSet<String>();
				m_go2proteinMap.put(complex, proteins);
			}
			proteins.add(protein);
		}
		debug_println("# complex database #proteins " + m_protein2goMap.size());
		debug_println("# complex database #complexs " + m_go2proteinMap.size());
		//debug_println("# mismatch_genes " + mismatch_genes.size());
		inputFile.close();

		double sum = 0;
		for (String tf : m_go2proteinMap.keySet()) {
			// debug_println("\t"+ tf + "\t"+ m_go2proteinMap.get(tf).size());
			sum += m_go2proteinMap.get(tf).size();
		}
		debug_println("# Average proteins per complex: "
				+ (sum / m_go2proteinMap.size()));

		return m_go2proteinMap;
	}
	
	static void enrich_mod_ppi(double p_th) throws IOException {
		List<String> lines = FileUtil.getFileAsList(m_mod_file);
		boolean start = false;
		int total = 0;
		int match = 0;
		for (String line: lines) {
			if (line.startsWith("1")) {
				start = true;
			}
			if (!start) {
				continue;
			}
			total ++;
			String[] items = line.split("\t");
		    String ppi_str = items[4];
		    String[] protein_str = ppi_str.split(", ");
		    Set<String> proteins = new HashSet<String>();
		    for (String p: protein_str) {
		    	proteins.add(p);
		    }
			StringBuffer output = new StringBuffer();
			Set<String> res = enrich(proteins, output, p_th, 5000, 2, "Complex");			
			debug_println("Proteins: "+ proteins);
			debug_println(""+ output);
			if (res.size() > 1) {
				match ++;
			}
		}
		debug_println("REP>\t"+ m_mod_file +"\t"+ p_th +"\t" + 
				total +"\t"+ match +"\t" + ((double)match/total));
		
	}

	static void enrich_mod_iii(double p_th) throws Exception {
		Map<String, String> known2protein = GeneAnnot.getKnown2protein(m_conn);
		List<String> lines = FileUtil.getFileAsList(m_mod_file);
		boolean start = false;
		int total = 0;
		int match = 0;
		for (String line: lines) {
			if (line.startsWith("1")) {
				start = true;
			}
			if (!start) {
				continue;
			}
			total ++;
			String[] items = line.split("\t");
		    String ppi_str = items[4];
		    String[] protein_str = ppi_str.split(", ");
		    Set<String> proteins = new HashSet<String>();
		    for (String p: protein_str) {
		    	proteins.add(known2protein.get(p));
		    }
			StringBuffer output = new StringBuffer();
			Set<String> res = enrich(proteins, output, p_th, 5000, 2, "Complex");			
			debug_println("Proteins: "+ proteins);
			debug_println(""+ output);
			if (res.size() > 1) {
				match ++;
			}
		}
		debug_println("REP>\t"+ m_mod_file +"\t"+ p_th +"\t" + 
				total +"\t"+ match +"\t" + ((double)match/total));
		
	}
	
}
