//此UDF的速度增量0.1，误差0.1内收敛
#include "udf.h"
#include "mem.h"
#include "dynamesh_tools.h"
#define a1 0
#define b1 0
#define turbine1 20
/*
#define a1 0
#define b1 0
#define a2 0.650
#define b2 1.125833
#define a3 0.650
#define b3 -1.125833
#define turbine1 20
#define turbine2 21
#define turbine3 22
*/
real T_max;
real T_min;
real cycle=1.0,time2=2,temp_omega=0,total_T=0;
int change = 1;
int lunci = 0;	//以用于区分各宏函数是否为第一次运行
FILE *fp1;			//建立一个文件操作指针
FILE *fp2;			//建立一个文件操作指针
FILE *fp3;			//建立一个文件操作指针

DEFINE_ZONE_MOTION(core1, omega, axis, origin, velocity, time, dtime)
{
	int n;
	real temp1,temp2,cal_T_average,cal_w;//为了方便进行扭矩-转速求解，设立中间变量，此时计算扭矩、转速均为正值；
	real optimal_T,optimal_w,relative_e,T_average=0,mytime = 0,w=0;
	real x[ND_ND];
	real NV_VEC(A); /*用来保存入口处每个面单元的面积*/
	real torque_z=0;
	Domain *domain; /*流体计算域*/
	Thread *threadsurface; /*surface指针*/
	face_t f;
	if(lunci==0) 
	{
		fp1=fopen("*omega1.txt","w+");
		fp2=fopen("T_average1.txt","w+");
		fp3=fopen("relative_e1.txt","w+");
	}
	mytime = time-time2;	//这是第n次循环所处的时间位置
	
	*omega = temp_omega;
	fprintf(fp1,"omega:%g\n",temp_omega);
	origin[0] = a1;
	origin[1] = b1;

	if(time < 2.0) //先计算2S，此时扭矩正常波动，方便后续计算
	{
		*omega = -5;
		mytime = 0;
	}
	else
	{
		mytime = time-time2;
	}
	Message("\n*mytime:%g\n",mytime);
	w = *omega;
	if(w<0)
	{
		cal_w=-w;//此处为了后面计算方便，改成正值，使用数值进行计算，不考虑方向
	}		
	else
	{
		cal_w=w;
	}
	Message("\ncal_w:%g\n",cal_w);
	n = (int)(0.5/dtime); //强制转化，给出累加次数,计算逻辑为：运行0.5秒钟，除以时间步，就可以得到步数
	Message("\n n:%d\n",n);
	
	if(change == 1)
	{
		T_max = 0;
		T_min = 1000;
		change = 0;
	}
	Message("\nchange2:%d\n",change);
	
	if(mytime > 0.5*cycle)
	{
		domain = Get_Domain(1); 
		threadsurface = Lookup_Thread(domain, turbine1);
		begin_f_loop(f,threadsurface) // 循环surface面单元
		{
			F_CENTROID(x,f,threadsurface);
			F_AREA(A, f, threadsurface);  // 计算单元的面积
			torque_z+= (F_P(f, threadsurface) * A[1])*(x[0]-a1)-(F_P(f, threadsurface) * A[0])*(x[1]-b1);/*此时扭矩的符号（方向）与转速一致*/
		}
		end_f_loop(f,threadsurface)
		Message("\ntorque_z:%g\n",torque_z); 
		if(mytime>0.5 && mytime < 1) total_T+= torque_z; //如果mytime处于我们关注的那个周期，那么就累加转矩
	}
	if(mytime >= 0.995*cycle && change==0)
	{
		Message("\ntotal_T:%g\n",total_T);
		T_average=total_T/n;
		fprintf(fp2,"T_average:%g\n",T_average);
		Message("\nT_average:%g\n",T_average);
		
		
		if (T_average<0)
		{
			cal_T_average=-T_average; //这里是为了方便以后计算相对误差等数值，所以强制转为正值
		}
		else
		{
			cal_T_average=T_average;
		}
		Message("\ncal_T_average:%g\n",cal_T_average);
		
		optimal_T = 0.01558*cal_w*cal_w+0.0002267*cal_w+0.0006039;//optimal_T为正值，实际为负
		Message("\noptimal_T:%g\n",optimal_T);
		
		temp2 = sqrt(cal_T_average);//temp2与optimal_T情况一致
		optimal_w=8.4*temp2;
		Message("\noptimal_w:%g\n",optimal_w);
		
		temp1 = 0.1*(optimal_w-cal_w)+cal_w;//计算下一阶段的转速
		Message("\ntemp1:%g\n",temp1);
		relative_e = (temp1-cal_w)/cal_w;
		if(relative_e<0)relative_e = -relative_e;
		Message("\nrelative_e:%g\n",relative_e);
		fprintf(fp3,"relative_e:%g\n",relative_e);
		if(relative_e>0.05)
		{
			Message("\ntemp1:%g\n",temp1);
			if(T_average<0)
			{
				*omega=-temp1;
			}
			else
			{
				*omega=-13;//当实际扭矩为正，说明风力机反转，此时强制风力机转速到最佳转速附近
			}
			//Message("\n*omega:%g\n",*omega);
		}
		change = 1;
		Message("\nchange2:%d\n",change);
		time2 = time;
		total_T = 0 ;
		Message("\ntime2:%g\n",time2);
	}
	temp_omega = *omega;
	lunci++;
}
