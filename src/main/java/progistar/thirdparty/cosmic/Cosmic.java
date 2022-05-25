package progistar.thirdparty.cosmic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;


class CosmicMutation {
	public String cosmicID = null;
	public String somaticInfo = null;
}

public class Cosmic {
	
	public Hashtable<String, ArrayList<CosmicMutation>> mutations = new Hashtable<String, ArrayList<CosmicMutation>>();
	
	public Cosmic (String[] paths, String[] formats) {
		try {
			for(int i = 0; i<paths.length; i++) {
				String path = paths[i];
				String format = formats[i];
				BufferedReader BR = new BufferedReader(new FileReader(new File(path)));
				String line = null;
				BR.readLine(); // skip header
				
				boolean isXPrint = false;
				boolean isYrint = false;
				boolean isMPrint = false;
				
				int cosmicIDIdx = 11;
				int somaticInfoIdx = 16;
				int gPosIdx = 27;
				
				if(format.equalsIgnoreCase("format2")) {
					cosmicIDIdx = 16;
					somaticInfoIdx = 30;
					gPosIdx = 38;
				}
				
				while((line = BR.readLine()) != null) {
					line = line.replaceAll("\t", ",");
					
					String[] fields = line.split(",");
					
					if(format.equalsIgnoreCase("format2")) {
						if(fields.length <= gPosIdx) {
							continue;
						}
					}
					
					String cosmicID = fields[cosmicIDIdx];
					String somaticInfo = fields[somaticInfoIdx];
					String gPos = fields[gPosIdx].replace("g.", "");
					gPos = gPos.replace("m.", "");
					
					if(!isXPrint && (gPos.contains("x") || gPos.contains("X"))) {
						System.out.println(gPos);
						isXPrint = true;
					}
					if(!isYrint && (gPos.contains("y") || gPos.contains("Y"))) {
						System.out.println(gPos);
						isYrint = true;
					}
					if(!isMPrint && (gPos.contains("m") || gPos.contains("M"))) {
						System.out.println(gPos);
						isMPrint = true;
					}
					
					gPos = "chr"+gPos;
					
					CosmicMutation cs = new CosmicMutation();
					cs.cosmicID = cosmicID;
					cs.somaticInfo = somaticInfo;
					
					ArrayList<CosmicMutation> csList = this.mutations.get(gPos);
					if(csList == null) {
						csList = new ArrayList<CosmicMutation>();
						mutations.put(gPos, csList);
					}
					
					csList.add(cs);
				}
				
				BR.close();
			}
		}catch(IOException ioe) {
			
		}
	}
	
	
	public static void main(String[] args) throws IOException {
		String fileName = "/Users/gistar/projects/pXg/COSMICMutations/CosmicNCV.tsv";
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		System.out.println(BR.readLine().replaceAll("\t", ","));
		System.out.println(BR.readLine().replaceAll("\t", ","));
		System.out.println(BR.readLine().replaceAll("\t", ","));
		System.out.println(BR.readLine().replaceAll("\t", ","));
		System.out.println(BR.readLine().replaceAll("\t", ","));
		System.out.println(BR.readLine().replaceAll("\t", ","));
		
		BR.close();
		
	}
}
