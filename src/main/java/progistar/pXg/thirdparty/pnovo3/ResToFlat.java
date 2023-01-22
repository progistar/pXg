package progistar.pXg.thirdparty.pnovo3;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

class PNOVORES {
	String scanStr = "";
	
	ArrayList<String> records = new ArrayList<String>();
	
}

public class ResToFlat {

	
	public static void main(String[] args) throws IOException {
		File pNovoResFile = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/pNovo3/S4.pNovo.res");
		BufferedReader BR = new BufferedReader(new FileReader(pNovoResFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(pNovoResFile.getAbsolutePath().replace(".res", ".flat.res")));
		
		String line = null;
		String header = "FileName\tRank\tPeptide\tScore";
		BW.append(header);
		BW.newLine();
		
		ArrayList<PNOVORES> resList = new ArrayList<PNOVORES>();
		PNOVORES res = null;
		while((line = BR.readLine()) != null) {
			if(line.startsWith("S")) {
				res = new PNOVORES();
				String scan = line.split("\t")[1];
				res.scanStr = scan;
				resList.add(res);
			} else if(line.startsWith("P")) {
				String[] field = line.split("\t");
				res.records.add(field[0].substring(1)+"\t"+field[1]+"\t"+field[2]);
			}
		}
		
		BR.close();
		
		
		for(PNOVORES res_ : resList) {
			for(String record : res_.records) {
				BW.append(res_.scanStr+"\t"+record);
				BW.newLine();
			}
		}
		
		BW.close();
		
	}
}
