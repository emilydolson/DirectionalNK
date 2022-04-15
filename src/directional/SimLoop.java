/**
 * Copyright (C) 2017, Sonia Singhal
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
package directional;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * @author Sonia
 *
 */
public class SimLoop {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int N = 15;
		double mu[] = new double [] {0.000005, (double)1/15};
		int Kvals[] = new int [] {0, 1, 2, 5, 10, 14};
		long seed;
		long sseed;
		// double cuts[] = new double [] {0.1, 10, 0.71, 100};
		int gens = 50;
		boolean getsShocks[] = new boolean [] {false, true};
		String [] argv = new String[28];
		String directory = "outs/";
		argv[0] = "-n";
		argv[1] = Integer.toString(N);
		argv[2] = "-k";
		argv[4] = "-s";
		argv[6] = "-o";
		argv[8] = "-l";
		argv[10] = "-c";
		argv[11] = "{0.4}";
		// System.out.println(argv[11]);
		argv[12] = "-r";
		//argv[13] = Double.toString(mu);
		argv[14] = "-g";
		argv[15] = Integer.toString(gens);
		argv[16] = "-f";
		argv[17] = "resources/population.txt";
		argv[18] = "-m";
		argv[19] = "multi_random";
		argv[20] = "-tracepop";
		argv[21] = "true";
		argv[22] = "-shocks";
		argv[24] = "-sseed";
		argv[26] = "-maxfit";
		argv[27] = "1.0";
		PrintStream out;
		try {
			out = new PrintStream(new File("files.txt"));
			//Run Sudden replicates, then Gradual replicates
			for(int sh = 0; sh < getsShocks.length; sh++) {
				argv[23] = (getsShocks[sh]) ? "{2 0.4000000 4 0.4095833 6 0.4191667 8 0.4287500 10 0.4383333 12 0.4479167 14 0.4575000 16 0.4670833 18 0.4766667 20 0.4862500 22 0.4958333 24 0.5054167 26 0.5150000 28 0.5245833 30 0.5341667 32 0.5437500 34 0.5533333 36 0.5629167 38 0.5725000 40 0.5820833 42 0.5916667 44 0.6012500 46 0.6108333 48 0.6204167 50 0.6300000}" : "{2 0.63 4 0.63 6 0.3 8 0.63 10 0.63 12 0.63 14 0.63 16 0.63 18 0.63 20 0.63 22 0.63 24 0.63 26 0.63 28 0.63 30 0.63 32 0.63 34 0.63 36 0.63 38 0.63 40 0.63 42 0.63 44 0.63 48 0.63 50 0.63 51 0.63}";
				//Run original mutation rate of 5x10^-6, then mutation rate of 1/genome
				for(int mut = 0; mut <mu.length; mut++) {
					argv[13] = Double.toString(mu[mut]);
					System.out.println(argv[13]);
			//Loop over replicates (20)
			for(int i = 0; i < 20; i++){
				//Loop over K values
				for(int k = 0; k < Kvals.length; k++){
					seed = (long)(Math.random()*50000);
					sseed = (long)(Math.random()*50000);
					String file = (getsShocks[sh]) ? directory+"N"+N+"-k"+Kvals[k]+"-s"+seed+"-r"+mu[mut]+"-Gradual.txt" : directory+"N"+N+"-k"+Kvals[k]+"-s"+seed+"-r"+mu[mut]+"-Sudden.txt";
					//String file = directory +N+"-"+k+"-"+seed+".txt";
					argv[3] = Integer.toString(Kvals[k]);
					argv[5] = Long.toString(seed);
					argv[7] = directory;
					argv[9] = file;
					argv[25] = Long.toString(sseed);
					write(out,file);
					// Configuration config = new Configuration(argv);
					// Simulation sim = new Simulation(argv);
					//sim.main(argv);
					Simulation.main(argv);
					System.out.println("done "+file);
				}

			}
			}
			}
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	static void write(PrintStream out,String line){
		out.println(line);
		return;
	}
}
