package algorithm;

import CloudSim.CloudSimExe;
import CloudSim.CloudSimPrint;
import Excel.WriteToExcel;
import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.core.CloudSim;

import java.text.DecimalFormat;
import java.util.*;

import static CloudSim.CloudSimExe.getCloudletById;
import static CloudSim.CloudSimPrint.getVmById;

public class MOEAISa {
    static List<Vm> vmList = CloudSimExe.vmList;
    static List<Cloudlet> cloudletList = CloudSimExe.cloudletList;
    static Random random = new Random(System.currentTimeMillis());
    private static int AGENT = 20;  // 种群大小
    private static int MAX_IT = 500;

    private static int[][] Population;
    private static double[] Fitness;
    private static int[] BestSol;
    private static double BestFitness;

    public static void runSimulationMOEAISa(DatacenterBroker broker) {
        initialize(cloudletList.size(), vmList.size());

        // 执行MOEAISa调度方案
        int[] schedule = MOEAISaScheduling();
        assignResourcesWithSchedule(schedule);

        CloudSim.startSimulation();
        List<Cloudlet> newList = broker.getCloudletReceivedList();
        CloudSim.stopSimulation();

        System.out.println("======================MOEAISa=======================");
        double finishTm = getMaxTimeOfSchedule(BestSol);
        System.out.println("最大完成时间：" + finishTm);

        double totalTime = getSumTimeOfSchedule(BestSol);
        System.out.println("系统执行时间：" + totalTime);

        System.out.println("资源利用率：" + totalTime / (finishTm * vmList.size()) * 100);
    }

    private static void initialize(int taskNum, int vmNum) {
        Population = new int[AGENT][taskNum];
        Fitness = new double[AGENT];
        BestSol = new int[taskNum];
        BestFitness = Double.MAX_VALUE;

        for (int i = 0; i < AGENT; i++) {
            for (int j = 0; j < taskNum; j++) {
                Population[i][j] = random.nextInt(vmNum);
            }
            Fitness[i] = fitness(Population[i]);
            if (Fitness[i] < BestFitness) {
                BestFitness = Fitness[i];
                System.arraycopy(Population[i], 0, BestSol, 0, taskNum);
            }
        }
    }

    public static int[] MOEAISaScheduling() {
        for (int iter = 0; iter < MAX_IT; iter++) {
            int[][] Parents = selectParents();
            int[][] Offspring = selfAdaptiveCrossover(Parents);
            mutation(Offspring);

            for (int i = 0; i < AGENT; i++) {
                double newFitness = fitness(Offspring[i]);
                if (newFitness < Fitness[i]) {
                    Fitness[i] = newFitness;
                    System.arraycopy(Offspring[i], 0, Population[i], 0, cloudletList.size());
                    if (newFitness < BestFitness) {
                        BestFitness = newFitness;
                        System.arraycopy(Offspring[i], 0, BestSol, 0, cloudletList.size());
                    }
                }
            }
        }
        return BestSol;
    }

    private static int[][] selectParents() {
        int[][] parents = new int[AGENT][cloudletList.size()];
        for (int i = 0; i < AGENT; i++) {
            int idx1 = random.nextInt(AGENT);
            int idx2 = random.nextInt(AGENT);
            parents[i] = Fitness[idx1] < Fitness[idx2] ? Population[idx1] : Population[idx2];
        }
        return parents;
    }

    private static int[][] selfAdaptiveCrossover(int[][] parents) {
        int[][] offspring = new int[AGENT][cloudletList.size()];
        for (int i = 0; i < AGENT; i++) {
            for (int j = 0; j < cloudletList.size(); j++) {
                double crossoverProb = random.nextDouble() < 0.5 ? 0.9 : 0.1;
                offspring[i][j] = random.nextDouble() < crossoverProb ? parents[i][j] : random.nextInt(vmList.size());
            }
        }
        return offspring;
    }

    private static void mutation(int[][] offspring) {
        double mutationRate = 0.1;
        for (int i = 0; i < AGENT; i++) {
            if (random.nextDouble() < mutationRate) {
                for (int j = 0; j < cloudletList.size(); j++) {
                    offspring[i][j] = random.nextInt(vmList.size());
                }
            }
        }
    }

    public static void assignResourcesWithSchedule(int[] schedule) {
        for (int i = 0; i < schedule.length; i++) {
            getCloudletById(i).setVmId(schedule[i]);
        }
    }

    private static double fitness(int[] schedule) {
        return getMaxTimeOfSchedule(schedule);
    }

    public static double getMaxTimeOfSchedule(int[] schedule) {
        double maxTime = 0;
        Map<Integer, ArrayList<Integer>> vmTasks = new HashMap<>();
        for (int i = 0; i < cloudletList.size(); i++) {
            vmTasks.computeIfAbsent((int) schedule[i], k -> new ArrayList<>()).add(i);
        }
        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }
            double runtime = length / getVmById(vmTask.getKey()).getMips();
            if (maxTime < runtime) {
                maxTime = runtime;
            }
        }
        return maxTime;
    }

    public static double getSumTimeOfSchedule(int[] schedule) {
        double sumTime = 0;
        Map<Integer, ArrayList<Integer>> vmTasks = new HashMap<>();
        for (int i = 0; i < cloudletList.size(); i++) {
            vmTasks.computeIfAbsent((int) schedule[i], k -> new ArrayList<>()).add(i);
        }
        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }
            double runtime = length / getVmById(vmTask.getKey()).getMips();
            sumTime += runtime;
        }
        return sumTime;
    }
}
