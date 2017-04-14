package edu.nchu.syslab.ppi;

import java.util.ArrayList;
import java.util.List;

import edu.usc.zhoulab.common.RocResult;

//import ccliu.umls.SVMResults;

import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;

/**
 * Runs an svm prediction on input nodes and return an SVMResults object
 * 
 */
public class SVMPredictor2 {

	public static RocResult svmPredict(svm_model model,
			List<svm_node[]> nodes, boolean expected) {
		RocResult results = new RocResult();
		for (svm_node[] n : nodes) {
			//double svm_predict = svm.svm_predict(model, n);
			double[] dec_values = new double[1];
			svm.svm_predict_values(model, n, dec_values);
			int label = (dec_values[0] > 0.0) ? 1 : -1;
			results.addSample(expected, label, dec_values[0]);
		}
		return results;
	}
	/*static double m_min_dec = 0;
	static double m_max_dec = 0;*/
	static List<Double> m_dec_list = null;
	public static RocResult svmPredict(svm_model model,
			List<svm_node[]> nodes, boolean expected, double threshold) {
		RocResult results = new RocResult();
		/*m_min_dec = 0;
		m_max_dec = 0;*/
		m_dec_list = new ArrayList<Double>();
		for (svm_node[] n : nodes) {
			//double svm_predict = svm.svm_predict(model, n);
			double[] dec_values = new double[1];
			svm.svm_predict_values(model, n, dec_values);
			//System.out.println("  svm_predict " + svm_predict+" -> "+ dec_values[0]);
			//int label = (svm_predict > 0.0) ? 1 : -1;
			int label = (dec_values[0] >= threshold) ? 1 : -1;
			m_dec_list.add(dec_values[0]);
			results.addSample(expected, label, dec_values[0]);
			/*if (dec_values[0] < m_min_dec){
				m_min_dec = dec_values[0];
			}
			if (dec_values[0] > m_max_dec){
				m_max_dec = dec_values[0];
			}*/
		}
		return results;
	}
	
	public static double svm_get_decval(svm_model model, svm_node[] x) {
		int nr_class = 2;

		double[] dec_values = new double[nr_class * (nr_class - 1) / 2];

		svm.svm_predict_values(model, x, dec_values);
		System.out.println("    svm_get_decval dec_values " + dec_values[0]);
		int pos2 = 0;
		double dv = 0;
		double decval = 0;
		double maxabsdecval = 0;
		for (int i = 0; i < nr_class; i++)
			for (int j = i + 1; j < nr_class; j++) {
				dv = dec_values[pos2++];
				if (Math.abs(dv) > maxabsdecval)
					maxabsdecval = Math.abs(dv);
				decval = dv;
			}

		return decval;
	}
	public static boolean svmPredict(svm_model model, svm_node[] n) {
		boolean label = (svm.svm_predict(model, n) > 0.0);
		return label;
	}
}
