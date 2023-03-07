package progistar.thirdparty.netMHCpan;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Hashtable;

import progistar.pXg.constants.Parameters;

public class NetMHCpanRunner {

	private NetMHCpanRunner() {};
	
	private static Hashtable<String, Boolean> peptides = new Hashtable<String, Boolean>();
	
	public static void putPeptide(String peptide) {
		peptides.put(peptide, true);
	}
	
	public static void run (String option) {
		try {
			StringBuilder command = new StringBuilder();
			
			command.append(Parameters.netMHCpanPath).append(" ").append(option);
			
			Process netMHCpan = Runtime.getRuntime().exec(command.toString());
			
			try(BufferedReader input = new BufferedReader(new InputStreamReader(netMHCpan.getInputStream()))) {
	            String line;

	            while ((line = input.readLine()) != null) {
	                System.out.println(line);
	            }
	        }

			
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String args[]) {
		run("");
	}
}
