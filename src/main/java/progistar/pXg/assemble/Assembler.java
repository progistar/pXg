package progistar.pXg.assemble;

import java.util.ArrayList;
import java.util.Hashtable;

import progistar.pXg.data.XBlock;

public class Assembler {
	private Assembler() {}

	private static ArrayList<XBlock> 			xBlocks		 	= new ArrayList<XBlock>();
	private static Hashtable<XBlock, String> 	xBlockToPSeq	= new Hashtable<XBlock, String>();
	
	/**
	 * Simply add XBlock to arraylist.<br>
	 * @param xBlock
	 */
	public static void addXBlock (String pSeq, XBlock xBlock) {
		xBlocks.add(xBlock);
		xBlockToPSeq.put(xBlock, pSeq);
	}
	
	public static String getPSeq (XBlock xBlock) {
		return xBlockToPSeq.get(xBlock);
	}
	
	public static ArrayList<XBlock> assemble () {
		ArrayList<XBlock> contigXBlocks = new ArrayList<XBlock>();
		
		// zero size list
		if(xBlocks.size() == 0) {
			return contigXBlocks;
		}
		
		Hashtable<String, ArrayList<XBlock>> hashXBlocks = new Hashtable<String, ArrayList<XBlock>>();
		for(XBlock xBlock : xBlocks) {
			// strand and genomic sequence are enough to represent singularity.
			String key = xBlock.strand+"_"+xBlock.genomicSequence;
			
			ArrayList<XBlock> thisXBlocks = hashXBlocks.get(key);
			if(thisXBlocks == null) {
				thisXBlocks = new ArrayList<XBlock>();
			}
			
			thisXBlocks.add(xBlock);
			hashXBlocks.put(key, thisXBlocks);
		}
		
		// get contigXBlocks
		hashXBlocks.forEach((key, xBlocks) -> {
			contigXBlocks.addAll(getContigs(xBlocks));
		});
		
		return contigXBlocks;
	}
	
	/**
	 * Assume that all Xblocks here have the same genomicseuqence for peptide with same strand. <br> 
	 * 
	 * @param xBlocks
	 * @return
	 */
	private static ArrayList<XBlock> getContigs (ArrayList<XBlock> xBlocks) {
		ArrayList<XBlock> contigXBlocks = new ArrayList<XBlock>();
		if(xBlocks.size() == 0) {
			return contigXBlocks;
		}
		
		ArrayList<ContigNode> contigNodes = new ArrayList<ContigNode>();
		
		// assume that all xBlocks here are same.
		for(int i=0; i<xBlocks.size(); i++) {
			XBlock xBlock = xBlocks.get(i);
			
			String fullSequence = "X"+xBlock.fullReadSequence+"X";
			String[] flankSequences = fullSequence.split(xBlock.genomicSequence);
			
			if(flankSequences.length > 2) {
				System.out.println("unexpected ... duplicated sequences at contig generator");
			}
			
			ContigNode contigNode = new ContigNode();
			contigNode.xBlocks.add(xBlock);
			contigNode.leftFlankSequence = flankSequences[0].replaceFirst("X", "");
			contigNode.rightFlankSequence = flankSequences[1].replaceFirst("X", "");
			contigNode.targetSequence = xBlock.genomicSequence;
			
			contigNodes.add(contigNode);
		}
		
		// merge contig nodes O(n2): do not care performance ... because the size of blocks is expected to be very small.
		int seed = 1;
		for(int i=0; i<contigNodes.size(); i++) {
			ContigNode seedContigNode = contigNodes.get(i);
			// not assigned yet
			if(seedContigNode.seed == -1) {
				// assign seed
				seedContigNode.seed = seed++;
				seedContigNode.isSeed = true;
				
				for(int j=i+1; j<contigNodes.size(); j++) {
					ContigNode contigNode = contigNodes.get(j);
					// this contig node already joins other family.
					if(contigNode.seed != -1) {
						continue;
					}

					// 1: seed is longer than it
					// -1: it is longer than seed
					// 0: different sequence
					byte left = 0;
					byte right = 0;
					
					// check left
					if(seedContigNode.leftFlankSequence.contains(contigNode.leftFlankSequence)) {
						left = 1;
					} else if(contigNode.leftFlankSequence.contains(seedContigNode.leftFlankSequence)) {
						left = -1;
					}
					
					// check right
					if(seedContigNode.rightFlankSequence.contains(contigNode.rightFlankSequence)) {
						right = 1;
					} else if(contigNode.rightFlankSequence.contains(seedContigNode.rightFlankSequence)) {
						right = -1;
					}
					
					// merge
					if(left != 0 && right != 0) {
						if(left == -1) {
							seedContigNode.leftFlankSequence = contigNode.leftFlankSequence;
						}
						if(right == -1) {
							seedContigNode.rightFlankSequence = contigNode.rightFlankSequence;
						}
						
						contigNode.seed = seedContigNode.seed;
						seedContigNode.xBlocks.addAll(contigNode.xBlocks);
					}
				}
			}
		}
		
		// only present xBlocks of seed contig nodes
		for(int i=0; i<contigNodes.size(); i++) {
			ContigNode contigNode = contigNodes.get(i);
			if(contigNode.isSeed) {
				String locus = contigNode.getLocus();
				String contig = contigNode.getFullSequence();
				
				// seed contig has at least one xBlock
				// and all xBlocks in the same contig must be the same.
				XBlock xBlock = contigNode.xBlocks.get(0);
				xBlock.fullReadSequence = null;
				xBlock.genomicLocus = locus;
				xBlock.genomicSequence = contig;
				xBlock.targetReadCount = contigNode.xBlocks.size();
				
				
				contigXBlocks.add(xBlock);
			}
		}
		
		return contigXBlocks;
	}
}
