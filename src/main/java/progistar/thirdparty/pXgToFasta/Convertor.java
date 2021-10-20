package progistar.thirdparty.pXgToFasta;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Convertor {

	public static void main(String[] args) throws IOException {
		String pXgResult = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\Laumont_NatCommun2016\\Results";
		BufferedReader BR = new BufferedReader(new FileReader(pXgResult));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
		}
		
		BR.close();
	}
}
