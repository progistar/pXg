package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class TargetDecoyStat {

	public static void main(String[] args) throws IOException {
		File file = new File("C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results\\2.withCalibration\\3.S3_pXg\\PeptideAnnotationS3_5ppm_002_recal.rep1.rank10.pXg.fdr.dist");
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		double cCutoffScore = 0;
		double ncCutoffScore = 0;
		
		BR.readLine(); // skip header
		double cTargets = 0;
		double cDecoys = 0;
		double ncTargets = 0;
		double ncDecoys = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			double score = Double.parseDouble(fields[0]);
			
			cTargets += Double.parseDouble(fields[1]);
			cDecoys += Double.parseDouble(fields[2]);
			ncTargets += Double.parseDouble(fields[3]);
			ncDecoys += Double.parseDouble(fields[4]);
			
			if(cTargets != 0) {
				if( cDecoys/cTargets < 0.05 ) {
					cCutoffScore = score;
				}
			}
			if(ncTargets != 0) {
				if( ncDecoys/ncTargets < 0.05 ) {
					System.out.println(score);
					ncCutoffScore = score;
				}
			}
		}
		
		System.out.println("c-threshold: "+cCutoffScore);
		System.out.println("nc-threshold: "+ncCutoffScore);
		
		BR.close();
	}
}
