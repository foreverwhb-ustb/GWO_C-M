package algorithm;

import CloudSim.CloudSimExe;
import CloudSim.CloudSimPrint;
import Excel.WriteToExcel;
import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.core.CloudSim;

import java.text.DecimalFormat;
import java.util.*;

import static CloudSim.CloudSimExe.getCloudletById;
import static CloudSim.CloudSimPrint.getVmById;

public class GTOHBA {

    // 算法参数
    private static int MAX_IT = 500;  // 最大迭代次数
    private static int AGENT = 20;     // 种群大小
    private static double p = 0.5;     // 概率参数
    private static double beta = 0.5;  // 调整参数

    static List<Vm> vmList = CloudSimExe.vmList;
    static List<Cloudlet> cloudletList = CloudSimExe.cloudletList;

    // 种群位置
    private static ArrayList<Double> population = new ArrayList<>();

    // 最佳个体位置和适应度
    private static double bestPosition;
    private static double bestFitness;

    private static Random random = new Random(System.currentTimeMillis());

    // 计算适应度（这里假设适应度函数是位置的平方，实际问题中需根据问题定义）
    private static double fitness(int[] schedule) {
        double fitness = 0;
        fitness = getMaxTimeOfSchedule(schedule);
        return fitness;
    }

    // 初始化种群
    private static void generateInitialPopulation() {
        for (int i = 0; i < AGENT; i++) {
            population.add(random.nextDouble());
        }
    }

    private static int[] getScheduleFromPosition(double position) {
        // 假设位置是一个浮动值，可以通过它来决定任务分配
        // 这里使用 position 来生成一个简单的调度方案
        // 你可以根据实际的种群编码方式修改这个实现

        // 创建一个与任务数量相同的调度数组
        int[] schedule = new int[CloudSimExe.cloudletList.size()];  // 假设任务数量等于 cloudletList 的大小

        // 根据位置生成调度，假设位置对应某个虚拟机的分配
        for (int i = 0; i < schedule.length; i++) {
            // 通过位置来决定任务分配的虚拟机
            // 这里使用位置作为虚拟机的索引（可以根据实际的编码方式进行调整）
            schedule[i] = (int) ((position * i) % CloudSimExe.vmList.size());  // 将位置映射到虚拟机编号
        }

        return schedule;
    }

    // 计算最佳个体
    private static void updateBest() {
        for (double pos : population) {
            int[] schedule = getScheduleFromPosition(pos);  // 假设这里从位置获取调度方案
            double fitnessValue = fitness(schedule);  // 使用适应度函数
            if (fitnessValue < bestFitness) {
                bestFitness = fitnessValue;
                bestPosition = pos;
            }
        }
    }

    // 计算C
    private static double calculateC(int t) {
        return beta * (1 - (double) t / MAX_IT);
    }

    // 计算L
    private static double calculateL(double C) {
        return C;
    }

    // 更新个体位置
    private static double updatePosition(double currentPosition, double alpha, double gamma, double C, double L) {
        double r1 = random.nextDouble();
        double r2 = random.nextDouble();
        double r3 = random.nextDouble();

        double GX = 0;
        if (r1 >= p) {
            GX = (1 - 0) * r2 + 0;
        } else if (r1 >= 0.5) {
            GX = (r3 - C) * currentPosition + L * 1;
        } else {
            GX = currentPosition - gamma * L * (currentPosition - bestPosition) + r3 * (currentPosition - bestPosition);
        }

        return GX;
    }

    // GTOHBA主算法
    public static void runGTOHBA() {
        generateInitialPopulation();
        bestFitness = Double.MAX_VALUE;  // 设置初始最优适应度为最大值
        updateBest();

        // 主循环
        for (int t = 0; t < MAX_IT; t++) {
            double alpha = 2 * Math.exp(-t / (double) MAX_IT);
            double gamma = Math.sin(2.5 - t / (double) MAX_IT);

            // 计算C和L
            double C = calculateC(t);
            double L = calculateL(C);

            // 更新种群中的每个个体
            for (int i = 0; i < AGENT; i++) {
                double currentPosition = population.get(i);
                double newPosition = updatePosition(currentPosition, alpha, gamma, C, L);
                population.set(i, newPosition);
            }

            // 更新最佳个体
            updateBest();
        }
    }


    // ========================== 调度指标计算 ==========================

    public static void calculateSchedulingMetrics(DatacenterBroker broker) {
        // 执行GTOHBA调度方案
        runGTOHBA();

        CloudSim.startSimulation();

        // 获取所有任务的调度结果
        List<Cloudlet> newList = broker.getCloudletReceivedList();

        CloudSim.stopSimulation();

        // 从 Cloudlet 列表中提取调度方案，构造 schedule 数组
        int[] schedule = new int[newList.size()];
        for (int i = 0; i < newList.size(); i++) {
            schedule[i] = newList.get(i).getVmId();  // 获取每个任务分配的虚拟机编号
        }

        // 计算最大完成时间、系统执行时间和资源利用率
        double maxCompletionTime = getMaxTimeOfSchedule(schedule);
        double totalExecutionTime = getSumTimeOfSchedule(schedule);

        // 打印调度指标
        System.out.println("========================GTOHBA Scheduling Metrics=============================");
        System.out.println("最大完成时间: " + maxCompletionTime);
        System.out.println("系统执行时间: " + totalExecutionTime);
        System.out.println("资源利用率: " + totalExecutionTime / (maxCompletionTime*vmList.size())*100);
    }

    public static double getMaxTimeOfSchedule(int[] schedule) {
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

    public static double getSumTimeOfSchedule(int[] schedule) {
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



    // 执行GTOHBA算法调度时的最大完成时间、系统执行时间和资源利用率的计算
    public static void runSimulationGTOHBA(DatacenterBroker broker) {
        calculateSchedulingMetrics(broker);
    }
}
