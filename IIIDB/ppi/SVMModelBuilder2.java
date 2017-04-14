package edu.nchu.syslab.ppi;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import edu.usc.zhoulab.common.RocResult;

//import ccliu.umls.SVMResults;

import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;
import libsvm.svm_parameter;
import libsvm.svm_problem;

/**
 * Build an SVM classifier with three options: LINEAR, RBF, PLOY.
 *
 */
public class SVMModelBuilder2
{
	//public static final int SVM_TYPE = svm_parameter.LINEAR;
	//public static final int SVM_TYPE = svm_parameter.POLY;
	//public static final int SVM_TYPE = svm_parameter.RBF;
	//public static int SVM_TYPE = 0;
	public static final int cMin = -10;
	public static final int cMax = 10;
	public static final int cStep = 2;
	public static final int gMin = -10;
	public static final int gMax = 10;
	public static final int gStep = 2;
	public static final int numCV = 4;
	static Random random = new Random();
	
	/*public static svm_model buildModel(List<svm_node[]> positives, List<svm_node[]> negatives)
	{
		// return buildModel(positives, negatives, SVM_TYPE, 1.0/(positives.size()+negatives.size()), 1);
		// default for manifold
		return buildModel(positives, negatives, SVM_TYPE, Math.pow(2,0), Math.pow(2,0));
	}*/
	
	public static svm_model buildModel(List<svm_node[]> positives, List<svm_node[]> negatives, int SVM_TYPE, double gamma, double C)
	{
		svm_problem problem = new svm_problem();
		
		int numP = positives.size();
		int numN = negatives.size();
		int numNodes = numP+numN;
		/*
		System.out.println("numP "+numP);
		System.out.println("numN "+numN);
		System.out.println("positives "+positives);
		System.out.println("positives.get(0) "+positives.get(0));
		System.out.println("positives.get(0).length "+positives.get(0).length); */
		problem.x = new svm_node[numNodes][positives.get(0).length];
		double[] y = new double[numNodes];
		for (int i=0; i<numP; i++)
		{
			problem.x[i]=positives.get(i);
			y[i]=1.0;
		}
		for (int i=0; i<numN; i++)
		{
			y[numP+i]=-1.0;
			problem.x[numP+i]=negatives.get(i);
		}
		problem.y = y;
		problem.l = numNodes;
		
		svm_parameter params = new svm_parameter();
		
		if (SVM_TYPE == svm_parameter.LINEAR)
		{ // setting for Manifold July 24, 2008
			params.svm_type=svm_parameter.C_SVC;
			params.kernel_type=svm_parameter.LINEAR;
			params.gamma=gamma;
			params.cache_size=100;
			params.eps=0.001;
			params.C = C;
			params.nr_weight=2;
			params.weight_label = new int[] {-1,1};
		}

		if (SVM_TYPE == svm_parameter.RBF)
		{ 
			params.svm_type=svm_parameter.C_SVC;
			params.kernel_type=svm_parameter.RBF;
			params.gamma=0.001; // 1/k, k is the number of attributes in the input data
			params.cache_size=100;
			params.eps=0.001;
			params.C = 1.0; // default 1
			params.nr_weight=2;
			params.weight_label = new int[] {-1,1};
		}
		
		if (SVM_TYPE == svm_parameter.POLY)
		{ // setting for haifeng's comment July 24, 2008			
			params.svm_type=svm_parameter.C_SVC;
			params.kernel_type= svm_parameter.POLY;
			params.degree=5;
			params.gamma=1.0/(positives.size()+negatives.size());
			params.cache_size=100;
			params.eps=0.001;
			params.C = 1;
			params.nr_weight=2;
			params.weight_label = new int[] {-1,1};
		} 
		
		//System.err.println("Number of positive samples: "+numP);
		//System.err.println("Number of negative samples: "+numN);
		double pWeight = (double)numN/(double)numP;
		//System.err.println("positive weight: "+pWeight);
		params.weight = new double[] {1.0,pWeight};
		//System.err.println(svm.svm_check_parameter(problem, params));
		return svm.svm_train(problem, params);
	}	
	/*public static svm_model buildModel_byCG2(List<svm_node[]> positives, List<svm_node[]> negatives, int bestC, int bestG)	
	{
		System.out.println("buildModel_byCG() bestC="+bestC+" bestG="+bestG);		
		svm_model finalModel= buildModel(positives, negatives, SVM_TYPE, Math.pow(2,bestG), Math.pow(2,bestC));
		RocResult results = new RocResult();
		results.add(SVMPredictor2.svmPredict(finalModel, positives, true));
		results.add(SVMPredictor2.svmPredict(finalModel, negatives, false));				
		System.out.println("  training SVMResults "+results);		
		return finalModel;
	}
	
	public static svm_model buildModel_byCG(List<List<svm_node[]>> positives, List<List<svm_node[]>> negatives, int bestC, int bestG)	
	{
		System.out.println("buildModel_byCG() bestC="+bestC+" bestG="+bestG);
		List<svm_node[]> pos = new ArrayList<svm_node[]>();
		List<svm_node[]> neg = new ArrayList<svm_node[]>();
		for (List<svm_node[]> l:positives)
			pos.addAll(l);
		for (List<svm_node[]> l:negatives)
			neg.addAll(l);
		
		svm_model finalModel= buildModel(pos, neg, SVM_TYPE, Math.pow(2,bestG), Math.pow(2,bestC));
		RocResult results = new RocResult();
		results.add(SVMPredictor2.svmPredict(finalModel, pos, true));
		results.add(SVMPredictor2.svmPredict(finalModel, neg, false));				
		System.out.println("  training SVMResults "+results);		
		return finalModel;
	}
	
	public static svm_model buildModel_LOOCV(List<List<svm_node[]>> positives, List<List<svm_node[]>> negatives)	
	{
		int bestC=0;
		int bestG=0;
		double bestAcc=0;
		for (int g=gMin; g<=gMax; g+=gStep)
		//for (int g=gMax; g>=gMin; g-=gStep)
		{
			if (bestAcc > 0.8) {
				break;
			}
			for (int c=cMin; c<=cMax; c+=cStep)
			//for (int c=cMax; c>=cMin; c-=cStep)			
			{
				System.out.println("Trying c="+c+" g="+g);
				RocResult results = new RocResult();
				for(int i=0; i<positives.size(); i++) {
					List<svm_node[]> t_pnodes = positives.get(i);
					List<svm_node[]> t_nnodes = negatives.get(i);
					List<svm_node[]> pnodes = new ArrayList<svm_node[]>();
					List<svm_node[]> nnodes = new ArrayList<svm_node[]>();
					for(int j=0; j<positives.size(); j++) {
						if (i==j){
							continue;
						}
						pnodes.addAll(positives.get(j));
						nnodes.addAll(negatives.get(j));
					}
					svm_model model = buildModel(pnodes, nnodes, SVM_TYPE, Math.pow(2,g), Math.pow(2,c));
					results.add(SVMPredictor2.svmPredict(model, t_pnodes, true));
					results.add(SVMPredictor2.svmPredict(model, t_nnodes, false));
					//System.out.println("  i="+i+" c="+c+" g="+g+" Acc="+results.getTrueAccuracy() );
				}
				System.out.println(results);
				double acc = results.getTrueAccuracy();
				//System.out.println("c="+c+" g="+g+" Acc="+acc );
				if (acc>bestAcc+0.0001)
				{
					bestC=c;
					bestG=g;
					bestAcc=acc;
				}				
			}
		}		
		System.out.println("bestC="+bestC+" bestG="+bestG+" bestAcc="+bestAcc );
		List<svm_node[]> pos = new ArrayList<svm_node[]>();
		List<svm_node[]> neg = new ArrayList<svm_node[]>();
		for (List<svm_node[]> l:positives)
			pos.addAll(l);
		for (List<svm_node[]> l:negatives)
			neg.addAll(l);
		
		svm_model finalModel= buildModel(pos, neg, SVM_TYPE, Math.pow(2,bestG), Math.pow(2,bestC));
		RocResult results = new RocResult();
		results.add(SVMPredictor2.svmPredict(finalModel, pos, true));
		results.add(SVMPredictor2.svmPredict(finalModel, neg, false));				
		System.out.println("  training SVMResults "+results);		
		return finalModel;
	}

	public static svm_model buildModel_LOOCV_byFDR(List<List<svm_node[]>> positives, List<List<svm_node[]>> negatives)	
	{
		int bestC=0;
		int bestG=0;
		double bestRate=1.1;
		for (int g=gMin; g<=gMax; g+=gStep)
		{
			if (bestRate < 0.2) {
				break;
			}
			for (int c=cMin; c<=cMax; c+=cStep)
			{
				System.out.println("Trying c="+c+" g="+g);
				RocResult results = new RocResult();
				for(int i=0; i<positives.size(); i++) {
					List<svm_node[]> t_pnodes = positives.get(i);
					List<svm_node[]> t_nnodes = negatives.get(i);
					List<svm_node[]> pnodes = new ArrayList<svm_node[]>();
					List<svm_node[]> nnodes = new ArrayList<svm_node[]>();
					for(int j=0; j<positives.size(); j++) {
						if (i==j){
							continue;
						}
						pnodes.addAll(positives.get(j));
						nnodes.addAll(negatives.get(j));
					}
					svm_model model = buildModel(pnodes, nnodes, SVM_TYPE, Math.pow(2,g), Math.pow(2,c));
					results.add(SVMPredictor2.svmPredict(model, t_pnodes, true));
					results.add(SVMPredictor2.svmPredict(model, t_nnodes, false));					
				}
				System.out.println(results);
				double rate = results.getFDR();
				System.out.println("  FDR="+rate);
				if (rate<bestRate-0.00001)
				{
					bestC=c;
					bestG=g;
					bestRate=rate;
				}				
			}
		}		
		System.out.println("bestC="+bestC+" bestG="+bestG+" bestRate="+bestRate );
		List<svm_node[]> pos = new ArrayList<svm_node[]>();
		List<svm_node[]> neg = new ArrayList<svm_node[]>();
		for (List<svm_node[]> l:positives)
			pos.addAll(l);
		for (List<svm_node[]> l:negatives)
			neg.addAll(l);
		
		svm_model finalModel= buildModel(pos, neg, SVM_TYPE, Math.pow(2,bestG), Math.pow(2,bestC));
		RocResult results = new RocResult();
		results.add(SVMPredictor2.svmPredict(finalModel, pos, true));
		results.add(SVMPredictor2.svmPredict(finalModel, neg, false));				
		System.out.println("  training SVMResults "+results);		
		return finalModel;
	}
	
	public static List<Integer> getRandomIndexes(List<Integer> indexes, int num)
	{
		List<Integer> ret = new ArrayList<Integer>();
		if (indexes.size()<=num)
			ret.addAll(indexes);
		else
		{
			while (num>0)
			{
				int index = random.nextInt(indexes.size());
				ret.add(indexes.get(index));
				indexes.set(index, indexes.get(indexes.size()-1));
				indexes.remove(indexes.size()-1);
				num--;
			}
		}
		return ret;
	}
	
	public static List<svm_node[]> getSubset(List<svm_node[]> nodes, int numSamples)
	{
		List<svm_node[]> ret = new ArrayList<svm_node[]>();
		if (nodes.size()<=numSamples)
			ret.addAll(nodes);
		else
		{
			int nodesSize=nodes.size();
			while (numSamples>0)
			{
				int i = random.nextInt(nodesSize);
				// add this to return array and exchange with last element
				svm_node[] myNode = nodes.get(i);
				ret.add(myNode);
				nodes.set(i, nodes.get(nodesSize-1));
				nodes.set(nodesSize-1, myNode);
				nodesSize--;
				numSamples--;
			}
		}
		return ret;
	}*/
		
}

