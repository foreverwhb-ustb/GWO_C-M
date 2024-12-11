package CloudSim;

import algorithm.*;
import algorithm.PSO.PSO;
import org.cloudbus.cloudsim.*;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.provisioners.BwProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.RamProvisionerSimple;

import java.io.*;
import java.io.File;
import java.util.*;

//import static algorithm.GA.runSimulationGA;

public class CloudSimExe {
    public static List<Cloudlet> cloudletList = new ArrayList<Cloudlet>();
    public static List<Vm> vmList;

    public static void RunTest(String dataFilePath, int taskNum) {

        long startTime = System.currentTimeMillis();

        Log.printLine("Starting to run simulations...");

        // 1. 创建实体之前初始化CloudSim工具包, 用户数量, 日历, 标志位
        int num_user = 1;
        Calendar calendar = Calendar.getInstance();
        boolean trace_flag = false;
        CloudSim.init(num_user, calendar, trace_flag);

        // 2. 创建数据中心, 数据中心是CloudSim的资源提供者.
        Datacenter datacenter0 = createDatacenter("Datacenter_0");
//        Datacenter datacenter1 = createDatacenter("Datacenter_1");

        // 3. 创建代理
        DatacenterBroker broker = createBroker();
        int brokerId = broker.getId();

        // 4. 创建虚拟机
        vmList = new ArrayList<Vm>();
        // 4.1 虚拟机描述
        int mips = 5000;
        long size = 100000; // image size (MB)
        int ram = 512; // vm memory (MB)
        long bw = 1000;
        int pesNumber = 1; // number of cpus
        String vmm = "Xen"; // VMM name

        // 4.2 创建五个虚拟机, 并添加到虚拟机队列

//        for (int i = 0; i < 10; i++) {
//            double tmpmips = 80 + (i % 3) * 10;
//            Vm vm = new Vm(i, brokerId, tmpmips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
//            vmList.add(vm);
//        }
        Vm vm1 = new Vm(0, brokerId, 500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm2 = new Vm(1, brokerId, 1000, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm3 = new Vm(2, brokerId, 1500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm4 = new Vm(3, brokerId, 500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm5 = new Vm(4, brokerId, 1000, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm6 = new Vm(5, brokerId, 1000, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm7 = new Vm(6, brokerId, 500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm8 = new Vm(7, brokerId, 1500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm9 = new Vm(8, brokerId, 500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm10 = new Vm(9, brokerId, 1000, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm11 = new Vm(10, brokerId, 500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        Vm vm12 = new Vm(11, brokerId, 1500, pesNumber, ram, bw, size, vmm, new CloudletSchedulerSpaceShared());
        vmList.add(vm1); vmList.add(vm2); vmList.add(vm3); vmList.add(vm4);
//        vmList.add(vm5); vmList.add(vm6); vmList.add(vm7); vmList.add(vm8);
//        vmList.add(vm9); vmList.add(vm10); vmList.add(vm11); vmList.add(vm12);

        // 4.3 将创建的虚拟机提交至代理
        broker.submitVmList(vmList);

        System.out.println(vmList.size());

        // 5. 创建任务提交至代理
        createTasks(brokerId, dataFilePath, taskNum);
        broker.submitCloudletList(cloudletList);

        System.out.println(cloudletList.size());

        // 6. 开始模拟
//        GWO.runSimulationGWO(broker);
//        MGWO.runSimulationGWO(broker);
//        main.java.algorithm.IBGWO.runSimulationGWO(broker);
//        TSMGWO.runSimulationGWO(broker);
        CMGWO2.runSimulationGWO(broker);
//        MOEAISa.runSimulationMOEAISa(broker);
//        SaPSO.runSimulationSaPSO(broker);
//        MOMRank.runSimulationMOMRank(broker);
//        MBO.runSimulationmBO(broker);
//        GTOHBA.runSimulationGTOHBA(broker);
        long endTime = System.currentTimeMillis();
        System.out.println("总时间：" + (endTime - startTime));

    }

    //    private static void createTasks(int brokerId, String filePath, int taskNum) {
//        try {
//            @SuppressWarnings("resource")
//            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
//            String data = null;
//            int index = 0;
//
////            File file = new File("data\\cloudlets.txt");//文件路径
////            FileReader fileReader = new FileReader(file);
////            LineNumberReader reader = new LineNumberReader(fileReader);
//
//            //cloudlet properties.
//            int pesNumber = 1;
//            long fileSize = 100;
//            long outputSize = 100;
//            UtilizationModel utilizationModel = new UtilizationModelFull();
//
//            while ((data = br.readLine()) != null) {
//                System.out.println(data);
//                String[] taskLength = data.split(" ");
//                for (int i = 0; i < taskLength.length; i++) {
//                    if (i == taskLength.length) {
//                        br.close();
//                        return;
//                    }
//                    Cloudlet task = new Cloudlet(i, (long) Double.parseDouble(taskLength[i]), pesNumber, fileSize, outputSize, utilizationModel, utilizationModel, utilizationModel);
//                    task.setUserId(brokerId);
//                    cloudletList.add(task);
//                }
//            }
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }
    private static void createTasks(int brokerId, String filePath, int taskNum) {
        try {
            @SuppressWarnings("resource")
//            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
//            String data = null;
//            int index = 0;
//
//            //cloudlet properties.
//            int pesNumber = 1;
//            long fileSize = 1000;
//            long outputSize = 1000;
//            UtilizationModel utilizationModel = new UtilizationModelFull();
//
//            while ((data = br.readLine()) != null) {
//                System.out.println("data==="+data);
//                String[] taskLength = data.split("\t");
////                String[] taskLength = data.split("\\s+");
//
//                for (int j = 0; j < 50; j++) {
////                for (int j = 0; j < taskLength.length; j++) {
//                    Cloudlet task = new Cloudlet(index + j, (long) Double.parseDouble(taskLength[j]), pesNumber, fileSize,
//                            outputSize, utilizationModel, utilizationModel,
//                            utilizationModel);
//                    task.setUserId(brokerId);
//                    cloudletList.add(task);
//                    if (cloudletList.size() == taskNum) {
//                        br.close();
//                        return;
//                    }
//                }
//                //20 cloudlets each line in the file cloudlets.txt.
//                index += 50;
//            }
            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
            String data = null;
            int index = 0;

            // cloudlet properties
            int pesNumber = 1;
            long fileSize = 1000;
            long outputSize = 1000;
            UtilizationModel utilizationModel = new UtilizationModelFull();

            while ((data = br.readLine()) != null) {
                // 如果读取到空行，跳过
                if (data.trim().isEmpty()) {
                    continue;
                }

                // 打印数据进行检查
                System.out.println("data===" + data);

                // 按制表符分割任务长度
//                String[] taskLength = data.split("\t");
                String[] taskLength = data.split("\\s");
                System.out.println("taskLength===" + taskLength.length);
                System.out.println("taskNum===" + taskNum);

                if (taskLength.length == 0) {
                    continue;  // 如果没有任务长度，跳过
                }

                // 创建任务
                for (int j = 0; j < taskLength.length; j++) {
                    try {
                        // 解析任务长度并创建任务
                        Cloudlet task = new Cloudlet(index + j, (long) Double.parseDouble(taskLength[j]), pesNumber, fileSize,
                                outputSize, utilizationModel, utilizationModel, utilizationModel);
                        task.setUserId(brokerId);
                        cloudletList.add(task);
                        if (cloudletList.size() == taskNum) {
                            br.close();
                            return;
                        }
                    } catch (NumberFormatException e) {
                        System.out.println("Invalid number format in task length: " + taskLength[j]);
                        continue;
                    }
                }
//                index += taskLength.length;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static DatacenterBroker createBroker() {
        DatacenterBroker broker = null;
        try {
            broker = new DatacenterBroker("Broker");
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        return broker;
    }

    private static Datacenter createDatacenter(String name) {
        List<Host> hostList = new ArrayList<Host>();
        List<Pe> peList = new ArrayList<Pe>();

        int mips = 5000;
        peList.add(new Pe(0, new PeProvisionerSimple(mips))); // need to store MIPS Rating

        mips = 2500;
        peList.add(new Pe(1, new PeProvisionerSimple(mips)));

        mips = 2500;
        peList.add(new Pe(2, new PeProvisionerSimple(mips)));

        mips = 1500;
        peList.add(new Pe(3, new PeProvisionerSimple(mips)));

        mips = 1000;
        peList.add(new Pe(4, new PeProvisionerSimple(mips)));

        int hostId = 0;
        int ram = 4096; // host memory (MB)
        long storage = 10000000; // host storage
        int bw = 10000;

        hostList.add(new Host(hostId, new RamProvisionerSimple(ram),
                new BwProvisionerSimple(bw), storage, peList,
                new VmSchedulerTimeShared(peList)));
        String arch = "x86"; // system architecture
        String os = "Linux"; // operating system
        String vmm = "Xen";
        double time_zone = 10.0; // time zone this resource located
        double cost = 3.0; // the cost of using processing in this resource
        double costPerMem = 0.05; // the cost of using memory in this resource
        double costPerStorage = 0.001; // the cost of using storage in this
        // resource
        double costPerBw = 0.001; // the cost of using bw in this resource

        //we are not adding SAN devices by now
        LinkedList<Storage> storageList = new LinkedList<Storage>();

        DatacenterCharacteristics characteristics = new DatacenterCharacteristics(
                arch, os, vmm, hostList, time_zone, cost, costPerMem,
                costPerStorage, costPerBw);

        // 6. Finally, we need to create a PowerDatacenter object.
        Datacenter datacenter = null;
        try {
            datacenter = new Datacenter(name, characteristics,
                    new VmAllocationPolicySimple(hostList), storageList, 0);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return datacenter;
    }

    public static Cloudlet getCloudletById(int id) {
        for (Cloudlet c : cloudletList) {
            if (c.getCloudletId() == id)
                return c;
        }
        return null;
    }

    public static void writeTxtAppend(String file, String conent) {
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file, true)));
            out.write(conent + "\r\n");
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
