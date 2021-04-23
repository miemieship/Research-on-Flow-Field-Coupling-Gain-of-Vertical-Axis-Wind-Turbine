#include "udf.h"
#include "mem.h"
#include "dynamesh_tools.h"
#define surfaceID 11  /*surface处的ID*/
#define a 0
#define b 0

#define cycle 1; //每次变换角速度的周期
//用于记录周期信息的全局变量
real k = 0;		//全局变量，用于记录循环轮次（角速度改变的次数）
real torqueZ = 0;	//全局变量，用于记录每个时间步的扭矩

//用于计算平均转矩的数组，弄成数组是因为，要存储不同转速下对应的平均转矩
real T_max_list[200];
real T_min_list[200];  //1000才对，先改成0
real T_average_list[200]; //记录平均转矩
real w_list[200]=0;		 //记录转速
real p_list[200]=0;		//记录风机功率
real w_max;				//我们需要寻找的最大转速

DEFINE_EXECUTE_AT_END(surface_force)
{
	int j;
	real x[ND_ND];	
	real NV_VEC(A); /*用来保存入口处每个面单元的面积*/
	real sum_T_A = 0.0; /*入口总的速度*/
	real sum_A = 0.0; /*入口总的面积*/
	real force_x=0;
	real force_y=0;
	real torque_x=0;
	real torque_y=0;
	real torque_z=0;
	real p1=0;
	real p2=0;
	real mass = 10;
	Domain *domain; /*流体计算域*/
	Thread *threadsurface; /*surface指针*/
	face_t f;
	//初始化数组
	if(k==0)
	{
		for(j=0;j<200;j++)
		{
			T_max_list[j]=0;
			T_min_list[j]=1000;
		}	
	}
	domain = Get_Domain(1); /*获得流体计算域，对于单相流，id值为1；对于多相流，id值为大于1的整数。多相流中的id值可以在Phase对话框中查看*/
	threadsurface = Lookup_Thread(domain, surfaceID); /*查找surface的指针*/
	begin_f_loop(f,threadsurface) // 循环surface面单元
	{
                F_CENTROID(x,f,threadsurface);
		F_AREA(A, f, threadsurface);  // 计算单元的面积
		force_x += F_P(f, threadsurface) * A[0];
                force_y += F_P(f, threadsurface) * A[1];
				p1=x[0];
                p2=x[1];
                torque_z+=(F_P(f, threadsurface) * A[0])*p2-(F_P(f, threadsurface) * A[1])*p1;
               	}
	end_f_loop(f,threadsurface)
		torqueZ = torque_z;
      	Message("\nwhen time step is %d,the force of surface is %g\n",force_x,torque_z);
}

DEFINE_ZONE_MOTION(core, omega, axis, origin, velocity, time, dtime)
{
	real T = 0;
	real mytime = 0;	
	mytime = time - cycle*k;	//这是第n次循环所处的时间位置
	T = torque_z;	//每次时间步都会更新扭矩

	if(mytime > cycle) k++;	//当mytime超过一个定转速周期时，就说明，进入到了下一轮次，因此对k++
	
	*omega = 1*k+1; /*角速度赋值，从1开始，今后每个定转速周期＋0.1rad/s*/
	w_list[k] = *omega;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	
	if(mytime > 0.5*cycle)
	{
		if(T_max_list[k]<T)
			T_max_list[k]=T;
		if(T_min_list[k]>T)
			T_min_list[k]=T;
		T_average_list[k]=(T_max_list[k] + T_min_list[k])/2;
		p_list[k] = T_average_list[k] * w_list[k];
		Message("\n\n转速：%g\n平均转矩：%g\n对应功率：%g\n\n",w_list[k],T_average_list[k],p_list[k]);
	}
}
	