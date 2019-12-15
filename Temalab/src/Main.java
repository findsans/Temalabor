import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.util.Iterator;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import sun.misc.Unsafe;




public class Main {

	public static void main(String[] args) {
		try {
			disableWarning();
			List<Kmer> kmers = new ArrayList<>();
			String bc = "";
			double m = 0;
			double sg = 0;
			
			File file = new File("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\simcalls_1569510406.579453.tsv_processed.xlsx");
			FileInputStream fis = new FileInputStream(file);
			XSSFWorkbook wb = new XSSFWorkbook(fis);
			XSSFSheet sheet = wb.getSheetAt(0);
			
			Iterator<Row> rowIt = sheet.iterator();
			Row row = rowIt.next();
			while(rowIt.hasNext()) {
				row = rowIt.next();

				Iterator<Cell> cellIt = row.cellIterator();
				
				Cell cell1 = cellIt.next();
				Cell cell2 = cellIt.next();
				Cell cell3 = cellIt.next();
				
				bc = cell1.getStringCellValue();
				m = cell2.getNumericCellValue();
				sg = cell3.getNumericCellValue();
				
				kmers.add(new Kmer(bc, sg, m));
			}
			wb.close();
			
			/*week 1 test*/
			System.out.println("initializing...");
			HMM hmm = new HMM(kmers);
			
			/*hmm.setUpTransitionMatrix();
			
			BufferedWriter writer = new BufferedWriter(new FileWriter("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\transitionMatrix.txt"));
			for(int i = 0; i < 4096; i++) {
				for(int j = 0; j < 4096; j++) {
					writer.write(hmm.getTransitionMatrix()[i][j] + "\n");
				}
			}
			writer.close();*/
			
			double[][] om = new double[4096][4096];
			double[][] tm = new double[4096][4096];
			InputStreamReader isr1 = new InputStreamReader(new FileInputStream("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\overlapMatrix.txt"));
			BufferedReader br1 = new BufferedReader(isr1);
			InputStreamReader isr2 = new InputStreamReader(new FileInputStream("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\transitionMatrix.txt"));
			BufferedReader br2 = new BufferedReader(isr2);
			for(int i = 0; i < 4096; i++) {
				for(int j = 0; j < 4096; j++) {
					om[i][j] = Double.parseDouble(br1.readLine());
					tm[i][j] = Double.parseDouble(br2.readLine());
				}
			}
			hmm.setOverlapMatrix(om);
			hmm.setTransitionMatrix(tm);
			br1.close();
			br2.close();
			System.out.println("finished");
			
			/*Kmer start = hmm.getKmerByBaseChain("AAACCC");
			Kmer kmer1 = hmm.getKmerByBaseChain("AACCCG");
			Kmer kmer2 = hmm.getKmerByBaseChain("AACCCT");
			double test = hmm.similarity(start, start, kmer1, kmer2);
			if(test < 0)
				System.out.println("No transition exsists between the given Kmers.");
			else{
				kmer1.print();
				kmer2.print();
				System.out.println("similarity: " + test);
			}*/
			
			/*week 2 test*/
			for(int i = 0; i < 1000; i++) {
				hmm.generateRandomPath();
				
				/*BufferedWriter writer = new BufferedWriter(new FileWriter("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\results.txt", true));
				writer.append(hmm.getStart() + "\t" + hmm.getEnd() + "\n");
				writer.close();*/
				//System.out.print("\nReference sequence of kmers(mu, sigma):\n");
				System.out.print((i+1) + ". mérés folyamatban...");
				//hmm.printPath(hmm.getRef());
				BufferedWriter writer = new BufferedWriter(new FileWriter("C:\\Users\\Vili\\Documents\\Egyetem\\Témalabor\\results.txt", true));
				writer.append("Result number " + (i+1) + ":\n");
				writer.close();
				hmm.printPathToFile(hmm.getRef(), 0, -1);
			
				//System.out.print("Most similar sequences found:\n\n");
				hmm.findSimilarPaths(hmm.getRef());
				System.out.println(" Kész.");
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void disableWarning() {
	    try {
	        Field theUnsafe = Unsafe.class.getDeclaredField("theUnsafe");
	        theUnsafe.setAccessible(true);
	        Unsafe u = (Unsafe) theUnsafe.get(null);

	        @SuppressWarnings("rawtypes")
			Class cls = Class.forName("jdk.internal.module.IllegalAccessLogger");
	        Field logger = cls.getDeclaredField("logger");
	        u.putObjectVolatile(cls, u.staticFieldOffset(logger), null);
	    } catch (Exception e) {
	        // ignore
	    }
	}

}
