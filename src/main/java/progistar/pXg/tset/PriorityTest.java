package progistar.pXg.tset;

import java.io.IOException;

import progistar.pXg.utils.Priority;

public class PriorityTest {
	
	public static void main(String[] args) throws IOException {
		
		System.out.println(Priority.getRegionPenalty("ENST00000557508.5(9-;sense;N;C)","-"));
		
		//ENST00000554812.5(6-3N;sense;N;C)
		System.out.println(Priority.getRegionPenalty("ENST00000554812.5(6-3N;sense;N;C)","-"));
	}
	
}
