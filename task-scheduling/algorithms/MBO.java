package algorithm;

import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.core.CloudSim;
import CloudSim.CloudSimExe;
import CloudSim.CloudSimPrint;

import java.util.*;

import static CloudSim.CloudSimExe.getCloudletById;
import static CloudSim.CloudSimPrint.getVmById;

public class MBO {

    // 种群大小
    private static int POPULATION_SIZE = 20;
    // 最大迭代次数
    private static int MAX_ITERATIONS = 500;
    // 搜索空间维度（任务数量或虚拟机数量等相关维度）
    private static int DIMENSION;
    static List<Vm> vmList = CloudSimExe.vmList;
    static List<Cloudlet> cloudletList = CloudSimExe.cloudletList;

    // Bonobo个体类
    private static class Bonobo {
        int[] position;
        double fitness;

        Bonobo(int dimension) {
            position = new int[dimension];
            fitness = 0;
        }
    }

    // 种群列表
    private static ArrayList<Bonobo> population = new ArrayList<>();
    // 最优Bonobo个体
    private static Bonobo bestBonobo;

    // 初始化种群
    private static void initializePopulation() {
        Random random = new Random();
        DIMENSION = cloudletList.size();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Bonobo bonobo = new Bonobo(DIMENSION);
            for (int j = 0; j < DIMENSION; j++) {
                // 随机初始化Bonobo个体的位置（任务分配方式，例如随机分配任务到虚拟机）
                bonobo.position[j] = random.nextInt(vmList.size());
            }
            // 计算初始适应度
            bonobo.fitness = calculateFitness(bonobo.position);
            population.add(bonobo);
        }
        // 初始化最优Bonobo个体
        bestBonobo = population.get(0);
        for (Bonobo bonobo : population) {
            if (bonobo.fitness < bestBonobo.fitness) {
                bestBonobo = bonobo;
            }
        }
    }

    // 更新种群
    private static void updatePopulation(int iteration) {
        // 计算转换因子TF
        double TF = Math.exp(-((double) iteration / MAX_ITERATIONS));
        // 计算控制随机化参数DR
        double DR = 2 * new Random().nextDouble() - 1;
        // 控制随机化方向CD（设定为1，可根据需要调整）
        int CD = 1;

        ArrayList<Bonobo> newPopulation = new ArrayList<>();
        for (Bonobo bonobo : population) {
            Bonobo newBonobo = new Bonobo(DIMENSION);
            // 选择交配策略（这里简单随机选择一种，实际可根据概率等更复杂方式选择）
            int matingStrategy = new Random().nextInt(4);
            if (matingStrategy == 0) {
                // 限制性交配
                newBonobo.position = restrictiveMating(bonobo, bestBonobo, TF, CD, DR);
            } else if (matingStrategy == 1) {
                // 滥交
                newBonobo.position = promiscuousMating(bonobo, bestBonobo, TF, CD, DR);
            } else if (matingStrategy == 2) {
                // 外族交配
                newBonobo.position = extraGroupMating(bonobo, TF, CD, DR);
            } else {
                // 配偶制交配
                newBonobo.position = consortshipMating(bonobo, bestBonobo, TF, CD, DR);
            }
            // 计算新个体适应度
            newBonobo.fitness = calculateFitness(newBonobo.position);
            newPopulation.add(newBonobo);
        }
        population = newPopulation;
        // 更新最优Bonobo个体
        for (Bonobo bonobo : population) {
            if (bonobo.fitness < bestBonobo.fitness) {
                bestBonobo = bonobo;
            }
        }
    }

    // 限制性交配策略（示例，需完善）
    private static int[] restrictiveMating(Bonobo bonobo, Bonobo bestBonobo, double TF, int CD, double DR) {
        int[] newPosition = new int[DIMENSION];
        Random random = new Random();
        for (int i = 0; i < DIMENSION; i++) {
            // 根据论文公式实现位置更新，注意整数运算规则
            newPosition[i] = bonobo.position[i] + (int) (TF * CD * DR * 1.4 * (bestBonobo.position[i] - bonobo.position[i]) +
                    TF * CD * DR * 1.45 * (-1) * (bonobo.position[i] - bonobo.position[random.nextInt(DIMENSION)]));
        }
        return newPosition;
    }

    // 滥交策略（示例，需完善）
    private static int[] promiscuousMating(Bonobo bonobo, Bonobo bestBonobo, double TF, int CD, double DR) {
        int[] newPosition = new int[DIMENSION];
        Random random = new Random();
        for (int i = 0; i < DIMENSION; i++) {
            // 根据论文公式实现位置更新，注意整数运算规则
            newPosition[i] = bonobo.position[i] + (int) (TF * CD * DR * 1.4 * (bestBonobo.position[i] - bonobo.position[i]) +
                    TF * CD * DR * 1.45 * (1) * (bonobo.position[i] - bonobo.position[random.nextInt(DIMENSION)]));
        }
        return newPosition;
    }

    // 外族交配策略（示例，需完善）
    private static int[] extraGroupMating(Bonobo bonobo, double TF, int CD, double DR) {
        int[] newPosition = new int[DIMENSION];
        Random random = new Random();
        for (int i = 0; i < DIMENSION; i++) {
            // 根据论文公式实现位置更新（这里简化示例，未考虑所有条件）
            newPosition[i] = bonobo.position[i] + (int) (TF * CD * DR * (Math.exp(random.nextDouble() * 2 + 2 * random.nextDouble() - 2 / random.nextDouble()) *
                    (vmList.size() - 1 - bonobo.position[i])));
        }
        return newPosition;
    }

    // 配偶制交配策略（示例，需完善）
    private static int[] consortshipMating(Bonobo bonobo, Bonobo bestBonobo, double TF, int CD, double DR) {
        int[] newPosition = new int[DIMENSION];
        Random random = new Random();
        for (int i = 0; i < DIMENSION; i++) {
            // 根据论文公式实现位置更新（这里简化示例，未考虑所有条件）
            if (random.nextDouble() <= 0.0035) {
                newPosition[i] = bonobo.position[i] + (int) ((1) * Math.exp(-random.nextDouble()) *
                        (bonobo.position[i] - bonobo.position[random.nextInt(DIMENSION)]));
            } else {
                newPosition[i] = bonobo.position[random.nextInt(DIMENSION)];
            }
        }
        return newPosition;
    }


    // 根据任务分配方案分配资源（类似GWO中的分配方法，需根据实际情况调整）
    private static void assignResourcesWithBestSchedule(int[] schedule) {
        for (int i = 0; i < schedule.length; i++) {
            cloudletList.get(i).setVmId(schedule[i]);
        }
    }

    // 执行mBO调度算法
    public static void runSimulationmBO(DatacenterBroker broker) {

        CloudSimPrint print = new CloudSimPrint();

        // 初始化种群
        initializePopulation();

        // 迭代更新
        for (int i = 0; i < MAX_ITERATIONS; i++) {
            updatePopulation(i);
        }

        // 根据最优Bonobo个体分配任务
        assignResourcesWithBestSchedule(bestBonobo.position);

        CloudSim.startSimulation();

        // 最终步骤：模拟结束后打印结果
        List<Cloudlet> newList = broker.getCloudletReceivedList();

        CloudSim.stopSimulation();

        // 打印相关结果（例如最大完成时间、系统执行时间、资源利用率等，可参考GWO代码中的打印方式）
        double finishTm = getMaxTimeOfSchedule(bestBonobo.position);
        System.out.println("最大完成时间：" + finishTm);

        double totalTime = getSumTimeOfSchedule(bestBonobo.position);
        System.out.println("系统执行时间：" + totalTime);

        System.out.println("资源利用率：" + totalTime / (finishTm * vmList.size()) * 100);
    }

    private static double calculateFitness(int[] schedule) {

        double fitness = 0;
        fitness = getMaxTimeOfSchedule(schedule);

        return fitness;
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
                i--;
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
}