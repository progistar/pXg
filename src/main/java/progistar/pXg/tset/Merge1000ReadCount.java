package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Merge1000ReadCount {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pXg/OriginData/S4.RAW.PEAKS.pxg.pval.dist");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		BufferedWriter BW = new BufferedWriter(new FileWriter(file.getAbsolutePath().replace(".dist", ".1000.dist")));
		
		String header = BR.readLine();
		BW.append(header);
		BW.newLine();
		String line = null;
		
		double[][] countExp = new double[16][1001];
		double[][] countMock = new double[16][1001];
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			int peptLen = Integer.parseInt(fields[0]);
			int curCount = Integer.parseInt(fields[1]);
			double curExp = Double.parseDouble(fields[2]);
			double curMock = Double.parseDouble(fields[3]);
			
			if(curCount >= 1000) {
				curCount = 1000;
			}
			countExp[peptLen][curCount] += curExp;
			countMock[peptLen][curCount] += curMock;
		}
		
		for(int i=8; i<=15; i++) {
			for(int j=0; j<=1000; j++) {
				BW.append(i+"\t"+j+"\t"+countExp[i][j]+"\t"+countMock[i][j]);
				BW.newLine();
			}
		}
		
		BR.close();
		BW.close();
	}
}
