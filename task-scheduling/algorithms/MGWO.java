package algorithm;

import CloudSim.CloudSimExe;
import CloudSim.CloudSimPrint;
import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.core.CloudSim;

import java.text.DecimalFormat;
import java.util.*;

import static CloudSim.CloudSimExe.getCloudletById;
import static CloudSim.CloudSimPrint.getVmById;

public class MGWO {

    /**
     * 1 建模：
     * 每只狼：一个可能的任务分配序列
     * 虚拟机：猎物
     * 更新任务分配：攻击猎物
     * 虚拟机负载高：被更多的狼攻击
     * 选择正确的虚拟机：计算距离并选择合适的猎物
     * 分配任务到低负载的虚拟机：寻找正确的猎物
     * 2 算法流程
     * 1 初始化，随机分配任务
     * 2 计算负载阈值，检测每台虚拟机的负载情况，underload 或者 overload
     * 3 将 overload 的虚拟机上的任务重新分配出去
     */


    static List<Vm> vmList = CloudSimExe.vmList;
    static List<Cloudlet> cloudletList = CloudSimExe.cloudletList;
    static Random random = new Random(System.currentTimeMillis());
    private static int AGENT = 20;      // 种群大小
    //    private static int DIMENSION = 2;
    private static int MAX_IT = 500;

    static int Bandwidth = 500;

    private static int alpha;
    private static int beta;
    private static int delta;

    private static ArrayList<int[]> wolves = new ArrayList<int[]>();    // 干活的狼群，wolves[i]表示一个任务分配方式，也就是一个解。
    private static int[] alphaPos = new int[AGENT];
    private static int[] betaPos = new int[AGENT];
    private static int[] deltaPos = new int[AGENT];
    private static final DecimalFormat df = new DecimalFormat("0.00");

    static double crossoverProb = 0.9;
    static double mutationRate = 0.2;

    /**
     * 计算负载阈值
     *
     * @param schedule
     * @return
     */
    public static double CalLoadThreshold(int[] schedule) {
        double loadThreshold = 0;

        double[] vmLoad = new double[vmList.size()];
        double sumLoad = 0;
        for (int i = 0; i < schedule.length; i++) {
            int vmIndex = schedule[i];
            vmLoad[vmIndex] += cloudletList.get(i).getCloudletLength() / vmList.get(vmIndex).getMips();//每个虚拟机的负载
            sumLoad += vmLoad[vmIndex];//总负载
        }

        double AL = sumLoad / vmList.size();
        double sum = 0;
        for (int i = 0; i < vmList.size(); i++) {
            sum += Math.pow((vmLoad[i] - AL), 2);
        }
        double standardDev = sum / vmList.size();
        loadThreshold = AL + standardDev;
        return loadThreshold;
    }


    public static void runSimulationGWO(DatacenterBroker broker) {

        CloudSimPrint print = new CloudSimPrint();

        //执行GWO调度方案


        CloudSim.startSimulation();

        applyGWOScheduling();


        // Final step: Print results when simulation is over
        List<Cloudlet> newList = broker.getCloudletReceivedList();
//        String finishTm = print.printCloudletList(newList);
        print.printCloudletList(newList);

        CloudSim.stopSimulation();

//        for (Vm vm : vmList) {
//            Log.printLine(String.format("vm id= %s ,mips = %s ", vm.getId(), vm.getMips()));
//        }
//        String finishTm = print.printCloudletList(newList);

        System.out.println("========================MGWO=============================");
        double finishTm = getMaxTimeOfSchedule(alphaPos);
        System.out.println("最大完成时间：" + finishTm);
//
        double totalTime = getSumTimeOfSchedule(alphaPos);
        System.out.println("系统执行时间：" + totalTime);
//
        System.out.println("资源利用率：" + totalTime / (finishTm * vmList.size()) * 100);

        //由于每次执行GWO调度算法的调度结果都不同(由于GWO过程中加入了随机性,甚至可能比RR还差),以下3行代码是取n次调度方案来计算GA结果的平均执行时间.
//        int n = 10;
//        double avgRuntime = getAvgRuntimeByGWOScheduling(n);
//        System.out.println(String.format("==============Printing the average running time GWO_GA schedule plans ===================\n" +
//                "Avg runtime of (n=%d) GWO_GA schedule plans:%.2f ms.", n, avgRuntime));

//        WriteToExcel.write(newList, 6);
    }

    private static double getAvgRuntimeByGWOScheduling(int times) {
        double sum = 0;
        for (int i = 0; i < times; i++) {
            int[] schedule = getScheduleByGWO();
            double tempTime = getMaxTimeOfSchedule(schedule);
            sum += tempTime;
            wolves.clear();
        }
        return sum / times;
    }

    public static void applyGWOScheduling() {
        int[] schedule = getScheduleByGWO();
        assignResourcesWithSchedule(schedule);
    }

    public static void assignResourcesWithSchedule(int[] schedule) {
        for (int i = 0; i < schedule.length; i++) {
            getCloudletById(i).setVmId(schedule[i]);
        }
    }

    private static int[] getScheduleByGWO() {
        initialize(cloudletList.size(), vmList.size());
//        schedules = GWO(schedules);
        return GWO();
    }

    private static int[] findBestSchedule(ArrayList<int[]> schedules) {
        double bestFitness = 1000000000;
        int bestIndex = 0;
        for (int i = 0; i < schedules.size(); i++) {
            int[] schedule = schedules.get(i);
            double fitness = fitness(schedule);
            if (bestFitness > fitness) {
                bestFitness = fitness;
                bestIndex = i;
            }
        }
        return schedules.get(bestIndex);
    }

    public static int[] GWO() {
        for (int i = 0; i < MAX_IT; i++) {
            update(i);
//            updateByGA();
            update_wolves();
        }
        return alphaPos;
    }

    private static void updateByGA() {
        for (int i = 0; i < AGENT; i++) {
            int[] parent_1, parent_2;
            int a = random.nextInt(3);
            if (a == 0) {
                parent_1 = betaPos;
                parent_2 = deltaPos;
            } else if (a == 1) {
                parent_1 = alphaPos;
                parent_2 = deltaPos;
            } else {
                parent_1 = alphaPos;
                parent_2 = betaPos;
            }

            // 交叉
            if (random.nextDouble() < crossoverProb) {
                int crossPosition = new Random().nextInt(cloudletList.size() - 1);
                for (int j = crossPosition + 1; j < cloudletList.size(); j++) {
                    int tmp = parent_1[j];
                    parent_1[j] = parent_2[j];
                    parent_2[j] = tmp;
                }
            }
            int[] tmp = fitness(parent_1) > fitness(parent_2) ? parent_2 : parent_1;
            wolves.set(i, fitness(wolves.get(i)) > fitness(tmp) ? tmp : wolves.get(i));

            // 变异
            tmp = wolves.get(i);
            if (random.nextDouble() < mutationRate) {

                int mutationIndex = new Random().nextInt(cloudletList.size());
                int newVmId = new Random().nextInt(vmList.size());
                while (tmp[mutationIndex] == newVmId) {
                    newVmId = new Random().nextInt(vmList.size());
                }
            }
            wolves.set(i, fitness(wolves.get(i)) > fitness(tmp) ? tmp : wolves.get(i));
        }
    }

    public static void update_wolves() {
        // 更新三只头狼
        double result;
        int aux_index_alpha;
        int aux_index_beta;

        for (int i = 0; i < AGENT; i++) {
            result = fitness(wolves.get(i));
            if (result < fitness(alphaPos)) {
                aux_index_alpha = alpha;
                aux_index_beta = beta;

                alphaPos = Arrays.copyOfRange(wolves.get(i), 0, cloudletList.size());
                betaPos = Arrays.copyOfRange(wolves.get(aux_index_alpha), 0, cloudletList.size());
                deltaPos = Arrays.copyOfRange(wolves.get(aux_index_beta), 0, cloudletList.size());

                alpha = i;
                beta = aux_index_alpha;
                delta = aux_index_beta;

            } else if ((result < fitness(betaPos)) && (result > fitness(alphaPos))) {
                aux_index_beta = beta;

                betaPos = Arrays.copyOfRange(wolves.get(i), 0, cloudletList.size());
                deltaPos = Arrays.copyOfRange(wolves.get(aux_index_beta), 0, cloudletList.size());

                beta = i;
                delta = aux_index_beta;
            } else if ((result < fitness(deltaPos)) && (result > fitness(betaPos))) {
                delta = i;
                deltaPos = Arrays.copyOfRange(wolves.get(i), 0, cloudletList.size());
            }
        }
    }

    private static void update(int iteration) {
        // 更新狼群
        double a = update_a(iteration);
        //System.out.println(a);
        //new java.util.Scanner(System.in).nextLine();
        for (int i = 0; i < AGENT; i++) {
            int[] tempWolf = new int[cloudletList.size()];
            for (int j = 0; j < cloudletList.size(); j++) {

                double r1 = random.nextDouble();
                double r2 = random.nextDouble();
                double a1 = (2 * a * r1) - a;
                double c1 = 2 * r2;
                double dAlpha = ((c1 * alphaPos[j]) - (wolves.get(i)[j]));
                double x1 = Math.abs(alphaPos[j] - (a1 * dAlpha));

                r1 = random.nextDouble();
                r2 = random.nextDouble();
                double a2 = (2 * a * r1) - a;
                double c2 = 2 * r2;
                double dBeta = ((c2 * betaPos[j]) - (wolves.get(i)[j]));
                double x2 = Math.abs(betaPos[j] - (a2 * dBeta));

                r1 = random.nextDouble();
                r2 = random.nextDouble();
                double a3 = (2 * a * r1) - a;
                double c3 = 2 * r2;
                double dDelta = ((c3 * deltaPos[j]) - (wolves.get(i)[j]));
                double x3 = Math.abs(deltaPos[j] - (a3 * dDelta));

                int newPos = (int) (x1 + x2 + x3) / 3;
                if (newPos < 0 || newPos >= vmList.size()) {
                    newPos = random.nextInt(vmList.size());
                }

                tempWolf[j] = newPos;
            }

            // 判断当前这只狼更新完毕后，适应度值是否更好
            wolves.set(i, fitness(wolves.get(i)) > fitness(tempWolf) ? tempWolf : wolves.get(i));
        }

    }

    private static void initialize(int taskNum, int vmNum) {
        // 将任务随机分配到虚拟机上
        for (int i = 0; i < AGENT; i++) {
            //data structure for saving a schedule：array,index of array are cloudlet id,content of array are vm id.
            int[] schedule = new int[taskNum];
            for (int j = 0; j < taskNum; j++) {
                schedule[j] = new Random().nextInt(vmNum);
            }
            wolves.add(schedule);
        }

        alpha = elitism(0);
        beta = elitism(fitness(wolves.get(alpha)));
        delta = elitism(fitness(wolves.get(beta)));

        alphaPos = wolves.get(alpha);
        betaPos = wolves.get(beta);
        deltaPos = wolves.get(delta);
    }

    private static int elitism(double found_minimum) {
        double minimum = fitness(wolves.get(0));
        int i, index_minimum = 0;
        for (i = 1; i < AGENT; i++) {
            double tempFitness = fitness(wolves.get(i));
            if (tempFitness < minimum && tempFitness > found_minimum) {
                index_minimum = i;
                minimum = tempFitness;
            }
        }
        return index_minimum;   //返回适应度值最优的结果的狼的序号给alpha
    }

    private static double get_random(double min, double max) {
        return random.nextDouble() * (max - min) + min;
    }

    private static double update_a(double iteration) {
        return 2 - 2 * ((Math.log(iteration / MAX_IT) - 1) / (Math.log(1) - 1));
    }

    private static double fitness(int[] schedule) {
        // 适应度值计算公式：fitness = w1 * loadBalance + w2 * normalizedRunTime


        double fitness = 0;
        fitness = getMaxTimeOfSchedule(schedule) + getMinCost(schedule);

        return fitness;
    }

    private static double getMinCost(int[] schedule) {
        double minCost = 0;
        HashMap<Integer, ArrayList<Integer>> vmTasks = new HashMap<Integer, ArrayList<Integer>>();
        int size = cloudletList.size();
        for (int i = 0; i < size; i++) {
            // 得到每个虚拟机分得的任务序号，存在 vmTasks Map中
            if (!vmTasks.keySet().contains(schedule[i])) {
                ArrayList<Integer> taskList = new ArrayList<Integer>();
                taskList.add(i);
                vmTasks.put(schedule[i], taskList);
            } else {
                vmTasks.get(schedule[i]).add(i);
            }
        }
        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }

            // 计算每个虚拟机完成分配的任务所用的时间, 任务总长度/虚拟机的Mips
            double cost = (double) length / Bandwidth;
            if (minCost > cost) {
                minCost = cost;
            }
        }
        return minCost;

    }

    public static double getSumTimeOfSchedule(int[] schedule) {
        // 得到某个调度方案的运行时间
        double sumTime = 0;
        HashMap<Integer, ArrayList<Integer>> vmTasks = new HashMap<Integer, ArrayList<Integer>>();
        int size = cloudletList.size();
        for (int i = 0; i < size; i++) {
            // 得到每个虚拟机分得的任务序号，存在 vmTasks Map中
            if (!vmTasks.keySet().contains(schedule[i])) {
                ArrayList<Integer> taskList = new ArrayList<Integer>();
                taskList.add(i);
                vmTasks.put(schedule[i], taskList);
            } else {
                vmTasks.get(schedule[i]).add(i);
            }
        }

        double utilizationSum = 0;
        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            double utilization = 0;
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }

            // 计算每个虚拟机完成分配的任务所用的时间, 任务总长度/虚拟机的Mips
            double runtime = length / getVmById(vmTask.getKey()).getMips();
            sumTime += runtime;
        }

        // 计算得到最后完成任务的虚拟机所用的时间并返回.
        return sumTime;
    }

    public static double getMaxTimeOfSchedule(int[] schedule) {
        // 得到某个调度方案的运行时间
        double maxTime = 0;
        HashMap<Integer, ArrayList<Integer>> vmTasks = new HashMap<Integer, ArrayList<Integer>>();
        int size = cloudletList.size();
        for (int i = 0; i < size; i++) {
            // 得到每个虚拟机分得的任务序号，存在 vmTasks Map中
            if (!vmTasks.keySet().contains(schedule[i])) {
                ArrayList<Integer> taskList = new ArrayList<Integer>();
                taskList.add(i);
                vmTasks.put(schedule[i], taskList);
            } else {
                vmTasks.get(schedule[i]).add(i);
            }
        }

        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }

            // 计算每个虚拟机完成分配的任务所用的时间, 任务总长度/虚拟机的Mips
            double runtime = length / getVmById(vmTask.getKey()).getMips();
            if (maxTime < runtime) {
                maxTime = runtime;
            }
        }

        // 计算得到最后完成任务的虚拟机所用的时间并返回.
        return maxTime;
    }

    public static double getMinTimeOfSchedule(int[] schedule) {
        // 得到某个调度方案的运行时间
        double minTime = 1000;
        HashMap<Integer, ArrayList<Integer>> vmTasks = new HashMap<Integer, ArrayList<Integer>>();
        int size = cloudletList.size();
        for (int i = 0; i < size; i++) {
            // 得到每个虚拟机分得的任务序号，存在 vmTasks Map中
            if (!vmTasks.keySet().contains(schedule[i])) {
                ArrayList<Integer> taskList = new ArrayList<Integer>();
                taskList.add(i);
                vmTasks.put(schedule[i], taskList);
            } else {
                vmTasks.get(schedule[i]).add(i);
            }
        }

        for (Map.Entry<Integer, ArrayList<Integer>> vmTask : vmTasks.entrySet()) {
            int length = 0;
            for (Integer taskId : vmTask.getValue()) {
                length += getCloudletById(taskId).getCloudletLength();
            }

            // 计算每个虚拟机完成分配的任务所用的时间, 任务总长度/虚拟机的Mips
            double runtime = length / getVmById(vmTask.getKey()).getMips();
            if (minTime > runtime) {
                minTime = runtime;
            }
        }

        // 计算得到最后完成任务的虚拟机所用的时间并返回.
        return minTime;
    }

}
