package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class TargetDecoyStat {

	public static void main(String[] args) throws IOException {
		
		File[] fileList = new File[1];
		fileList[0] = new File("C:\\Users\\progi\\Desktop\\Projects\\pXg\\hello.M.out");
		
		BufferedReader BR = new BufferedReader(new FileReader(fileList[0]));
		String line = null;
		BR.readLine();
		
		Hashtable<String, ArrayList<String[]>> records = new Hashtable<String, ArrayList<String[]>> ();
		
		ArrayList<String> keys = new ArrayList<String>();
		while((line = BR.readLine()) != null) {
			String[] field = line.split("\t");
			String key = field[0]+"_"+field[1]+"_"+field[4];
			
			ArrayList<String[]> fields = records.get(key);
			if(fields == null) {
				fields = new ArrayList<String[]>();
				keys.add(key);
				records.put(key, fields);
			}
			fields.add(field);
		}
		BR.close();
		
		
		int targetIndex = 23;
		int decoyIndex = 24;
		int scoreIndex = 7;
		
		for(String key : keys) {
			ArrayList<String[]> fields = records.get(key);
			String[] best = null;
			for(String[] field : fields) {
				int tCount = Integer.parseInt(field[targetIndex]);
				int dCount = Integer.parseInt(field[decoyIndex]);
				
				if(tCount >= 6) {
					
					best = field;
					
					break;
				}
			}
			
			if(best != null) {
				for(int i=0; i<best.length; i++) {
					if(i!=0) System.out.print("\t");
					System.out.print(best[i]);
				}
				System.out.println();
			}
		}
	}
}
