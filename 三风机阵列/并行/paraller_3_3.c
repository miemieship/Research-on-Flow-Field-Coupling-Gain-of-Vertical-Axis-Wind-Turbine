#include "udf.h"
#include "mem.h"
#include "dynamesh_tools.h"

#define a3 0.650
#define b3 -1.125833

#define turbine3 22

real T_max3;
real T_min3;
real cycle3=1.0,time23=2,temp_omega3=0,total_T3=0;
int change3 = 1;
int lunci3 = 0;	//以用于区分各宏函数是否为第一次运行
real time_list3[100];
real omega_list3[100];
real T_average_list3[100];
real relative_e_list3[100];

DEFINE_ZONE_MOTION(core3, omega, axis, origin, velocity, time, dtime)
{
	int n,j;
	real temp1,temp2,cal_T_average,cal_w;//为了方便进行扭矩-转速求解，设立中间变量，此时计算扭矩、转速均为正值；
	real optimal_T,optimal_w,relative_e,T_average=0,mytime=0,w=0;
	real x[ND_ND];
	real NV_VEC(A); /*用来保存入口处每个面单元的面积*/
	real torque_z=0;
	real torque_z_sum = 0;
	Domain *domain; /*流体计算域*/
	Thread *threadsurface; /*surface指针*/
	face_t f;
	if(lunci3==0 && change3==1) 
	{
		for(j=0;j<100;j++)
		{
			time_list3[j]=0;
			omega_list3[j]=0;
			T_average_list3[j]=0;
			relative_e_list3[j]=0;
		}
	}
	mytime = time-time23;	//这是第n次循环所处的时间位置
	*omega = temp_omega3;
	time_list3[lunci3] = time;
	omega_list3[lunci3] = temp_omega3;
	origin[0] = a3;
	origin[1] = b3;

	if(time < 2.0) //先计算2S，此时扭矩正常波动，方便后续计算
	{
		*omega = -5;
		mytime = 0;
	}
	else
	{
		mytime = time-time23;
	}
	w = *omega;
	if(w<0)
	{
		cal_w = -w;//此处为了后面计算方便，改成正值，使用数值进行计算，不考虑方向
	}		
	else
	{
		cal_w = w;
	}
	n = (int)((cycle3-0.5)/dtime); //强制转化，给出累加次数,计算逻辑为：运行0.5秒钟，除以时间步，就可以得到步数
	
	if(change3 == 1)
	{
		T_max3 = 0;
		T_min3 = 1000;
		change3 = 0;
	}
	if(mytime > 0.5*cycle3)
	{
		domain = Get_Domain(1); 
		threadsurface = Lookup_Thread(domain, turbine3);
		begin_f_loop(f,threadsurface) // 循环surface面单元
		{
			F_CENTROID(x,f,threadsurface);
			F_AREA(A, f, threadsurface);  // 计算单元的面积
			torque_z+= (F_P(f, threadsurface) * A[1])*(x[0]-a3)-(F_P(f, threadsurface) * A[0])*(x[1]-b3);/*此时扭矩的符号（方向）与转速一致*/
		}
		end_f_loop(f,threadsurface)
		torque_z_sum = PRF_GRSUM1(torque_z);
		if(mytime>0.5 && mytime < 1) total_T3+= torque_z_sum; //如果mytime处于我们关注的那个周期，那么就累加转矩
	}

	if(mytime >= 0.995*cycle3 && change3==0)
	{
		T_average=total_T3/n;
		node_to_host_real_1(T_average);
		#if !RP_NODE
		T_average_list3[lunci3]=T_average;
		#endif
		
		if (T_average<0)
		{
			cal_T_average=-T_average; //这里是为了方便以后计算相对误差等数值，所以强制转为正值
		}
		else
		{
			cal_T_average=T_average;
		}
		
		optimal_T = 0.01073*cal_w*cal_w+0.481*cal_w-1.128;//optimal_T为正值，实际为负
		temp2 = cal_T_average + (optimal_T - cal_T_average)*0.3; //0.3为松弛因子,temp2为cal_T_average与optimal_T的中值，后面用TEMP2求出的optimal_w更加稳定
		optimal_w = -0.03706*(temp2*temp2)+1.784*temp2+2.338;//使用temp2来求转速
		temp1 = 0.3*(optimal_w-cal_w)+cal_w;//计算下一阶段的转速，0.3为松弛因子
		
		relative_e = (temp1-cal_w)/cal_w;
		
		if(relative_e<0)relative_e = -relative_e;
		relative_e_list3[lunci3]=relative_e;
		if(relative_e>0.01)
		{
			if(T_average<0)
			{
				*omega=-temp1;
			}
			else
			{
				*omega=-13;//当实际扭矩为正，说明风力机反转，此时强制风力机转速到最佳转速附近
			}
		}
		else
		{
			cycle3 = 2.0;
		}
		change3 = 1;
		time23 = time;
		total_T3 = 0 ;
		lunci3++;
	}
	temp_omega3 = *omega;
#if !RP_NODE
	Message("turbine3\n");
	for(j=0;j<lunci3;j++)
	{
		Message("time:%g\t\tomega:%g\t\tT_average:%g\t\trelative_e:%g\n",time_list3[j],omega_list3[j],T_average_list3[j],relative_e_list3[j]);
	}
#endif
}

