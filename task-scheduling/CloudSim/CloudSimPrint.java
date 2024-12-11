package CloudSim;

import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Vm;

import java.text.DecimalFormat;
import java.util.List;

public class CloudSimPrint {
    static List<Vm> vmList = CloudSimExe.vmList;
    public static String printCloudletList(List<Cloudlet> list) {
        int size = list.size();
        Cloudlet cloudlet;

        String indent = "    ";
        Log.printLine();
        Log.printLine("================ Execution Result ==================");
        Log.printLine("No." + indent + "Cloudlet ID" + indent + "STATUS" + indent
                + "Data center ID" + indent + "VM ID" + indent + "VM mips" + indent + "CloudletLength" + indent + "Time"
                + indent + "Start Time" + indent + "Finish Time");

        DecimalFormat dft = new DecimalFormat("###.##");
        for (int i = 0; i < size; i++) {
            cloudlet = list.get(i);
            Log.print(i + 1 + indent + indent + cloudlet.getCloudletId() + indent + indent);


            if (cloudlet.getStatus() == Cloudlet.SUCCESS) {
                Log.print("SUCCESS");

                Log.printLine(indent + indent + indent + cloudlet.getResourceId()
                        + indent + indent + indent + cloudlet.getVmId()
                        + indent + indent + getVmById(cloudlet.getVmId()).getMips()
                        + indent + indent + cloudlet.getCloudletLength()
                        + indent + indent + indent + indent
                        + dft.format(cloudlet.getActualCPUTime()) + indent
                        + indent + dft.format(cloudlet.getExecStartTime())
                        + indent + indent
                        + dft.format(cloudlet.getFinishTime()));
            }
        }
        Log.printLine("================ Execution Result Ends here ==================");
        //最后完成的任务的完成时刻就是调度方案的总执行时间
        return dft.format(list.get(size - 1).getFinishTime());
    }
    public static Vm getVmById(int vmId) {
        for (Vm v : vmList) {
            if (v.getId() == vmId)
                return v;
        }
        System.out.println("Warning: Vm with ID " + vmId + " not found in vmList.");
        return null;
    }
}
