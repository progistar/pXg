package progistar.pXg.assemble;

import java.util.ArrayList;

import progistar.pXg.data.XBlock;

public class ContigNode {

	public String leftFlankSequence		= null;
	public String targetSequence		= null;
	public String rightFlankSequence	= null;
	public ArrayList<XBlock> xBlocks	= null;
	public int seed						= -1;
	boolean isSeed						= false;
	
	public ContigNode () {
		this.xBlocks = new ArrayList<XBlock>();
	}
	
	public String getLocus () {
		int left = 0;
		int right = 0;
		
		String fullSequence = getFullSequence();
		
		left = fullSequence.indexOf(targetSequence) + 1;
		right = left + targetSequence.length() - 1;
		
		return "*:"+left+"-"+right;
	}
	
	public String getFullSequence () {
		StringBuilder fullSequence = new StringBuilder();
		
		fullSequence.append(leftFlankSequence).append(targetSequence).append(rightFlankSequence);
		
		return fullSequence.toString();
	}
}
