package algorithm;

import CloudSim.CloudSimExe;
import CloudSim.CloudSimPrint;
import Excel.WriteToExcel;
import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.core.CloudSim;

import java.util.*;

import static CloudSim.CloudSimExe.getCloudletById;
import static CloudSim.CloudSimPrint.getVmById;

public class SaPSO {
    static List<Vm> vmList = CloudSimExe.vmList;
    static List<Cloudlet> cloudletList = CloudSimExe.cloudletList;
    static Random random = new Random(System.currentTimeMillis());

    private static final int AGENT = 20;       // Population size
    private static final int MAX_IT = 500;     // Max iterations
    private static final double Vmax = 6.0;    // Max velocity
    private static final double wMax = 0.9;    // Max inertia weight
    private static final double wMin = 0.2;    // Min inertia weight
    private static final double c1 = 2;        // Cognitive parameter
    private static final double c2 = 2;        // Social parameter
    private static final double alpha = 0.1;   // Strategy probability update factor
    private static final double pMax = 0.9;    // Max strategy probability
    private static final double pMin = 0.1;    // Min strategy probability

    private static double[] strategyProb = {1 / 3.0, 1 / 3.0, 1 / 3.0}; // Strategy probabilities
    private static double[][] Population, Velocity;
    private static double[] pBestScore, gBestScore;
    private static double[] gBest;
    private static double[] cg_curve;

    public static void runSimulationSaPSO(DatacenterBroker broker) {
        initialize(cloudletList.size(), vmList.size());

        int[] schedule = SaPSOScheduling();
        assignResourcesWithSchedule(schedule);

        CloudSim.startSimulation();
        List<Cloudlet> newList = broker.getCloudletReceivedList();
        CloudSim.stopSimulation();

        System.out.println("======================SaPSO=======================");
        double finishTm = getMaxTimeOfSchedule(gBest);
        System.out.println("最大完成时间：" + finishTm);

        double totalTime = getSumTimeOfSchedule(gBest);
        System.out.println("系统执行时间：" + totalTime);

        System.out.println("资源利用率：" + totalTime / (finishTm * vmList.size()) * 100);
    }

    public static void assignResourcesWithSchedule(int[] schedule) {
        for (int i = 0; i < schedule.length; i++) {
            int vmId = schedule[i] % vmList.size(); // 强制限定ID范围
            Vm vm = getVmById(vmId);
            if (vm != null) {
                getCloudletById(i).setVmId(vmId);
            } else {
                System.out.println("Warning: Vm with ID " + vmId + " not found in vmList.");
                getCloudletById(i).setVmId(random.nextInt(vmList.size())); // 使用有效ID
            }
        }
    }

    private static void initialize(int taskNum, int vmNum) {
        Population = new double[AGENT][taskNum];
        Velocity = new double[AGENT][taskNum];
        pBestScore = new double[AGENT];
        gBest = new double[taskNum];
        gBestScore = new double[taskNum];
        cg_curve = new double[MAX_IT];

        Arrays.fill(pBestScore, Double.MAX_VALUE);
        Arrays.fill(gBestScore, Double.MAX_VALUE);

        for (int i = 0; i < AGENT; i++) {
            for (int j = 0; j < taskNum; j++) {
                Population[i][j] = random.nextInt(vmNum);
                Velocity[i][j] = random.nextDouble() * 2 * Vmax - Vmax; // Initialize velocity
            }
        }
    }

    public static int[] SaPSOScheduling() {
        for (int iter = 0; iter < MAX_IT; iter++) {
            for (int i = 0; i < AGENT; i++) {
                int[] schedule = new int[Population[i].length];
                for (int j = 0; j < schedule.length; j++) {
                    schedule[j] = (int) Population[i][j];
                }

                double fitness = fitness(Population[i]);
                if (pBestScore[i] > fitness) {
                    pBestScore[i] = fitness;
                    System.arraycopy(Population[i], 0, gBest, 0, Population[i].length);
                }
                if (gBestScore[i] > fitness) {
                    gBestScore[i] = fitness;
                    System.arraycopy(Population[i], 0, gBest, 0, Population[i].length);
                }
            }

            double w = wMax - iter * ((wMax - wMin) / MAX_IT);

            for (int i = 0; i < AGENT; i++) {
                String strategy = rouletteWheelSelection(strategyProb);
                switch (strategy) {
                    case "standard":
                        for (int j = 0; j < Population[i].length; j++) {
                            Velocity[i][j] = w * Velocity[i][j] + c1 * random.nextDouble() * (gBest[j] - Population[i][j])
                                    + c2 * random.nextDouble() * (gBest[j] - Population[i][j]);
                            Velocity[i][j] = Math.min(Math.max(Velocity[i][j], -Vmax), Vmax);
                            Population[i][j] += Velocity[i][j];
                        }
                        break;
                    case "mutation":
                        for (int j = 0; j < Population[i].length; j++) {
                            Population[i][j] += 0.1 * random.nextGaussian();
                        }
                        break;
                    case "crossover":
                        for (int j = 0; j < Population[i].length; j++) {
                            Population[i][j] = random.nextDouble() < 0.5 ? gBest[j] : Population[i][j];
                        }
                        break;
                }
            }

            cg_curve[iter] = Arrays.stream(gBestScore).min().orElse(Double.MAX_VALUE);
            updateStrategyProbabilities();
        }

        int[] finalSchedule = new int[gBest.length];
        for (int i = 0; i < finalSchedule.length; i++) {
            finalSchedule[i] = (int) gBest[i];
        }
        return finalSchedule;
    }

    private static void updateStrategyProbabilities() {
        double improvement = Math.max(0, random.nextDouble() * (1 - Arrays.stream(gBestScore).min().orElse(Double.MAX_VALUE)));
        for (int i = 0; i < strategyProb.length; i++) {
            strategyProb[i] = Math.min(Math.max(strategyProb[i] + alpha * improvement, pMin), pMax);
        }
        double sum = Arrays.stream(strategyProb).sum();
        for (int i = 0; i < strategyProb.length; i++) {
            strategyProb[i] /= sum; // Normalize
        }
    }

    private static String rouletteWheelSelection(double[] prob) {
        double r = random.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < prob.length; i++) {
            cumulative += prob[i];
            if (r <= cumulative) {
//                return switch (i) {
//                    case 0 -> "standard";
//                    case 1 -> "mutation";
//                    default -> "crossover";
//                };
                return "mutation";
            }
        }
        return "standard";
    }

    private static double fitness(double[] schedule) {
        return getMaxTimeOfSchedule(schedule);
    }

    public static double getMaxTimeOfSchedule(double[] schedule) {
        double maxTime = 0;
        Map<Integer, ArrayList<Integer>> vmTasks = new HashMap<>();

        for (int i = 0; i < cloudletList.size(); i++) {
            int vmId = (int) schedule[i];

            // 强制确保 vmId 在范围内
            vmId = Math.min(Math.max(vmId, 0), vmList.size() - 1);

            Vm vm = getVmById(vmId);
            if (vm == null) {
                System.out.println("Warning: Vm with ID " + vmId + " not found in vmList.");
                continue; // 跳过不存在的虚拟机
            }

            vmTasks.computeIfAbsent(vmId, k -> new ArrayList<>()).add(i);
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

    public static double getSumTimeOfSchedule(double[] schedule) {
        double sumTime = 0;
        Map<Integer, ArrayList<Integer>> vmTasks = new HashMap<>();
        for (int i = 0; i < cloudletList.size(); i++) {
            int vmId = (int) schedule[i] % vmList.size(); // 强制限制 vmId 在有效范围内
            Vm vm = getVmById(vmId);

            if (vm == null) {
                System.out.println("Warning: Vm with ID " + vmId + " not found in vmList.");
                continue; // 如果 vm 为 null，跳过此任务
            }

            vmTasks.computeIfAbsent(vmId, k -> new ArrayList<>()).add(i);
        }

        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }

            Vm vm = getVmById(vmTask.getKey());
            if (vm != null) {
                double runtime = length / vm.getMips();
                sumTime += runtime;
            }
        }
        return sumTime;
    }
}
