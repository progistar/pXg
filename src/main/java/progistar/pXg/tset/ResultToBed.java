package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

public class ResultToBed {

	public static void main(String[] args) throws IOException {
		String resultFileName = "";
		
		BufferedReader BR = new BufferedReader(new FileReader(resultFileName));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			
		}
		
		BR.close();
	}
}
