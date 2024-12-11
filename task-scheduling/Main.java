import CloudSim.CloudSimExe;
import CloudSim.CloudSimExe1;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Random;

public class Main {
    public static void main(String[] args) throws IOException {

//        System.out.println("Hello world!");
//        String testData = "data\\task5.txt";
        String testData = "data\\task4.txt";
//        createTestData(testData);
//        createTestData2();
        int taskNum = 200;
        CloudSimExe.RunTest(testData, taskNum);

//        File file = new File("data\\cloudlets.txt");//文件路径
//        FileReader fileReader = new FileReader(file);
//        LineNumberReader reader = new LineNumberReader(fileReader);
//        int number = 2;//设置指定行数
//        String txt = "";
//        int lines = 0;
//        while (txt != null) {
//            lines++;
//            txt = reader.readLine();
//            if (lines == number) {
//                System.out.println("第" + reader.getLineNumber() + "的内容是：" + txt + "\n");
//                System.exit(0);
//            }
//        }
//        reader.close();
//        fileReader.close();

    }

   public static void createTestData(String filePath) {
       //create 50,000 data as cloudlet length for subsequent testing.
       int taskNum = 300;
       int[] taskLength = new int[taskNum];
       for (int i = 0; i < taskNum; i++) {
//            taskLength[i] = (new Random().nextInt(200) + 1) * 50 + new Random().nextInt(500);
           taskLength[i] = new Random().nextInt(500, 1500);
       }

       StringBuilder sb = new StringBuilder();
       for (int i = 0; i < taskNum; i++) {
           sb.append(taskLength[i]).append("\t");
           if (i % 50 == 49)//20 data each line.
           {
               CloudSimExe1.writeTxtAppend(filePath, sb.toString());
               sb = null;
               sb = new StringBuilder();
           }
       }
   }
   public static void createTestData2() {
       int i = 500;
       int index = 1;
       String filePath = "data\\cloudlets.txt";
       StringBuilder sb = new StringBuilder();
       while (i <= 515) {
           int j = i;
           while (j < 1400) {
               sb.append(j).append("\t");
               j += 18;
               if (index % 50 == 0) {
                   CloudSimExe1.writeTxtAppend(filePath, sb.toString());
                   sb = null;
                   sb = new StringBuilder();
//                    System.out.println();
               }
               index++;
           }
           i += 3;
       }
   }
}
