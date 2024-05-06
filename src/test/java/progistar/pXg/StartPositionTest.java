package progistar.pXg;

import java.util.ArrayList;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import junit.framework.Assert;
import progistar.pXg.constants.Constants;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.TBlock;
import progistar.pXg.data.parser.GTFParser;

public class StartPositionTest {
	@DisplayName("Unmapped reads에 대한 테스트가 필요해")
	@Test
	public void test() {
		GenomicAnnotation annotation = GTFParser.parseGTF("/Users/gistar/resources/Models/gencode.v42.basic.annotation.gtf");


		// + strand for start
		TBlock tBlock = annotation.getTBlockByTID("ENST00000218388.9");
		System.out.println("dist from start: "+distFromStartSite(47583413, 47583439, tBlock, true));
		System.out.println("dist from stop: "+distFromStopSite(47583413, 47583439, tBlock, true));
		Assert.assertTrue(distFromStartSite(47583413, 47583439, tBlock, true).equalsIgnoreCase("-3~+24"));

		tBlock = annotation.getTBlockByTID("ENST00000223095.5");
		System.out.println("dist from start: "+distFromStartSite(101128394, 101128420, tBlock, true));
		System.out.println("dist from stop: "+distFromStopSite(101128394, 101128420, tBlock, true));
		Assert.assertTrue(distFromStartSite(101128394, 101128420, tBlock, true).equalsIgnoreCase("+1~+27"));
		// - strand for start
		tBlock = annotation.getTBlockByTID("ENST00000452404.6");
		System.out.println("dist from start: "+distFromStartSite(196941942, 196941971, tBlock, false));
		System.out.println("dist from stop: "+distFromStopSite(196941942, 196941971, tBlock, false));
		Assert.assertTrue(distFromStartSite(196941942, 196941971, tBlock, false).equalsIgnoreCase("-13~+17"));

		// + strand for stop
		tBlock = annotation.getTBlockByTID("ENST00000692248.1");
		System.out.println("dist from start: "+distFromStartSite(32696617, 32696643, tBlock, true));
		System.out.println("dist from stop: "+distFromStopSite(32696617, 32696643, tBlock, true));
		Assert.assertTrue(distFromStopSite(32696617, 32696643, tBlock, true).equalsIgnoreCase("+11~+37"));

		tBlock = annotation.getTBlockByTID("ENST00000545580.1");
		System.out.println("dist from start: "+distFromStartSite(60850293, 60850322, tBlock, true));
		System.out.println("dist from stop: "+distFromStopSite(60850293, 60850322, tBlock, true));
		Assert.assertTrue(distFromStopSite(60850293, 60850322, tBlock, true).equalsIgnoreCase("-30~-1"));

		// - strand for stop
		tBlock = annotation.getTBlockByTID("ENST00000614790.5");
		System.out.println("dist from start: "+distFromStartSite(38802084, 38802107, tBlock, false));
		System.out.println("dist from stop: "+distFromStopSite(38802084, 38802107, tBlock, false));
		Assert.assertTrue(distFromStopSite(38802084, 38802107, tBlock, false).equalsIgnoreCase("-13~+11"));

	}

	private String distFromStartSite (int start, int end, TBlock tBlock, boolean strand) {
		ArrayList<Byte> targets = new ArrayList<>();
		targets.add(Constants.CDS); targets.add(Constants.UTR5); targets.add(Constants.UTR3);
		int peptStart = tBlock.getRelativeLengthOfPosition(start, targets);
		int peptEnd = tBlock.getRelativeLengthOfPosition(end, targets);
		int startSite = tBlock.getStartSite();

		// return "-"
		// when 1) intron, 2) noncoding, 3) intergenic, 4) antisense, 5) unmapped
		if(peptStart == -1 || peptEnd == -1 || startSite == -1 || (strand != tBlock.strand)) {
			return "-";
		}

		if(strand) {
			peptStart = peptStart - startSite;
			peptEnd = peptEnd - startSite;
		} else {
			peptStart = startSite - peptStart;
			peptEnd = startSite - peptEnd;
			int swap = peptStart;
			peptStart = peptEnd;
			peptEnd = swap;
		}


		if(peptStart >= 0) {
			peptStart++;
		}
		if(peptEnd >= 0) {
			peptEnd++;
		}

		String out = "";
		if(peptStart > 0) {
			out += "+";
		}
		out += peptStart +"~";
		if(peptEnd > 0) {
			out += "+";
		}
		out += peptEnd;

		return out;
	}

	private String distFromStopSite (int start, int end, TBlock tBlock, boolean strand) {
		ArrayList<Byte> targets = new ArrayList<>();
		targets.add(Constants.CDS); targets.add(Constants.UTR5); targets.add(Constants.UTR3);
		int peptStart = tBlock.getRelativeLengthOfPosition(start, targets);
		int peptEnd = tBlock.getRelativeLengthOfPosition(end, targets);
		int endSite = tBlock.getStopSite();

		// return "-"
		// when 1) intron, 2) noncoding, 3) intergenic, 4) antisense, 5) unmapped
		if(peptStart == -1 || peptEnd == -1 || endSite == -1 || (strand != tBlock.strand)) {
			return "-";
		}

		if(strand) {
			peptStart = peptStart - endSite;
			peptEnd = peptEnd - endSite;
		} else {
			peptStart = endSite - peptStart;
			peptEnd = endSite - peptEnd;
			int swap = peptStart;
			peptStart = peptEnd;
			peptEnd = swap;
		}

		if(peptStart <= 0) {
			peptStart--;
		}
		if(peptEnd <= 0) {
			peptEnd--;
		}

		String out = "";
		if(peptStart > 0) {
			out += "+";
		}
		out += peptStart +"~";
		if(peptEnd > 0) {
			out += "+";
		}
		out += peptEnd;

		return out;
	}
}
