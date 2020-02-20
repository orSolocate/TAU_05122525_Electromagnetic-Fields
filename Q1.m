clear all;
clc;
close all;
    
simulate(0,1,(3),1)


  function simulate(yL,yH,SIGMA_1,q_num)
    h=11; % nodes
    l=1.0; % length
    w=1.0; % width
    x=linspace(0,l,h);
    y=linspace(0,w,h);
    V=zeros(h);
    sigma_0= 3; %the conductivity
    Sigma_1=SIGMA_1;
    YL=yL;
    YH=yH;
    YH_idx= find(((abs(y-YH))==min(abs((y-YH)))),1,'first');
    YL_idx= find(((abs(y- YL))==min(abs((y- YL)))),1,'first');
    N=5+1; % N<9 so : [xL, xH] = [0.2, 0.4].
    XL=0.2;
    XH=0.4;
    XH_idx= find((abs((x-XH))==min(abs((x-XH)))));
    XL_idx= find((abs((x- XL))==min(abs((x- XL)))));

    %%
    %iterations:
    for sigma_1 = Sigma_1
        V=zeros(h);
        % Boundary conditions:
        V(YL_idx:YH_idx,1)=0; % left
        V(YL_idx:YH_idx,h)=1; % right
        eps=1e-6;
        error=1;
        k=0;
        while error>eps && k<1000
            k=k+1;
            Vold=V;
            for i=2:h-1;
                for j=2:h-1;
                V(i,j)=.25*(V(i+1,j)+V(i-1,j)+V(i,j-1)+V(i,j+1));
                end
            end
            error=max(max(abs(Vold-V)));

            %Neumann Boundary conditions:
            V(h,1:h)=V(h-1,1:h); % Top 
            V(1,1:h)=V(2,1:h); % Bottom
            %conditions for the rest of the right and left sides - under and above the electrode: 
            V(YH_idx+1:h,1)=V(YH_idx+1:h,2); % above the left electrode.
            V(YH_idx+1:h,h)=V(YH_idx+1:h,h-1); % above the right electrode.
            V(1:YL_idx-1,1)=V(1:YL_idx-1,2); % under the left electrode.
            V(1:YL_idx-1,h)=V(1:YL_idx-1,h-1); % under the right electrode.
            if Sigma_1 ~= (3)
                %%conditions due to difference in conductivity:    
                V(1:h,XL_idx)=(1/(sigma_0+sigma_1))*(sigma_0*V(1:h,XL_idx-1)+(sigma_1*V(1:h,XL_idx+1))); %left of sigma_1
                V(1:h,XH_idx)=(1/(sigma_0+sigma_1))*(sigma_0*V(1:h,XH_idx+1)+(sigma_1*V(1:h,XH_idx-1))); %right of sigma_1    
            end
        end

        %Electric Field:
        [Ex,Ey]=gradient(V); %[C/V*m]
        Ex = -h*Ex; %we mult by h in order to correct the units. now:[C/V]
        Ey = -h*Ey;
        E =sqrt(Ex.^2+Ey.^2);

        %% figures as in 1c:
        figure('Renderer', 'painters', 'Position', [10 10 900 1600],'Name', sprintf('question %.0f', q_num))
        if Sigma_1 ~= (3)
            sgtitle(['sigma1 = ',num2str(sigma_1),'[S/m]'])
        end
        subplot(4,2,1),contourf(x,y,V),colormap;
        xlim([0, 1]),  ylim([0, 1])
        hold on;
        if Sigma_1 ~= (3)
            rectangle('Position',[XL 0 (XH-XL) 1],'EdgeColor','r');
        end
        title('Potential lines(steady-state)'),xlabel('x[m]'),ylabel('y[m]'),colorbar;
        hold off;

        subplot(4,2,2),pcolor(x,y,V),shading interp,colormap;
        xlim([0, 1]),  ylim([0, 1])
        hold on;
        if Sigma_1 ~= (3)
            rectangle('Position',[XL 0 (XH-XL) 1],'EdgeColor','r');
        end
        title('potential map (steady-state)'),xlabel('x[m]'),ylabel('y[m]'),colorbar;
        hold off;

        subplot(4,2,3),quiver(x,y,Ex,Ey,2);
        xlim([0, 1]),  ylim([0, 1])
        hold on;
        if Sigma_1 ~= (3)
            rectangle('Position',[XL 0 (XH-XL) 1],'EdgeColor','r');
        end
        title('E Field (steady-state)'),xlabel('x[m]'),ylabel('y[m]')
        hold off;

        subplot(4,2,4),contour(x,y,E); colorbar;
        xlim([0, 1]),  ylim([0, 1])
        hold on;
        if Sigma_1 ~= (3)
            rectangle('Position',[XL 0 (XH-XL) 1],'EdgeColor','r');
        end
        title('E Field Magnititude (steady-state)'),xlabel('x[m]'),ylabel('y[m]')
        hold off;

        y_05=find((abs((x-0.5))==min(abs((x-0.5)))),1,'first');
        y_03=find((abs((x-0.3))==min(abs((x-0.3)))),1,'first');

        subplot(4,2,5), plot(x,Ex(y_05,:));
        hold on;
        if Sigma_1 ~= (3)
            x_rect=x(XL_idx:XH_idx);
            plot(x_rect,Ex(y_05,XL_idx:XH_idx),'red')
        end
        title('Ex for y=0.5 (steady-state)'),xlabel('x[m]'),ylabel('Ex(y=0.5) [v/m]')
        hold off;

        subplot(4,2,6), plot(x,Ey(y_05,:));
        hold on;
        if Sigma_1 ~= (3)
            plot(x_rect,Ey(y_05,XL_idx:XH_idx),'red')
        end
        title('Ey for y=0.5 (steady-state)'),xlabel('x[m]'),ylabel('Ey(y=0.5) [v/m]')
        hold off;

        subplot(4,2,7), plot(x,Ex(y_03,:));
        hold on;
        if Sigma_1 ~= (3)
            plot(x_rect,Ex(y_03,XL_idx:XH_idx),'red')
        end
        title('Ex for y=0.3 (steady-state)'),xlabel('x[m]'),ylabel('Ex(y=0.3) [v/m]')
        hold off;

        subplot(4,2,8), plot(x,Ey(y_03,:));
        hold on;
        if Sigma_1 ~= (3)
            plot(x_rect,Ey(y_03,XL_idx:XH_idx),'red')
        end
        title('Ey for y=0.3 (steady-state)'),xlabel('x[m]'),ylabel('Ey(y=0.3) [v/m]')
        hold off;

    %%  
        %capacitance:

        % first, we shall calculate the charge Q on the left electrode (equals to the charge of the right electrode with opposite sign)
        % we will calculate the charge by using the boundery condition on the
        % perpendicular electric field:
        perpendicular_E1 = mean(Ex(YL_idx:YH_idx,h));
        perpendicular_E2 = 0; % assuming the fields outside the body is negligible
        charge_density= (perpendicular_E2- perpendicular_E1); % normalized by epsilon_zero.
        Q = charge_density*(YH-YL)*1; %1 is the unit length in the z direction.
        delta_V = 1; %the potential difference between the elctrodes.
        C = Q/ delta_V;

    %%
        % conductivity:

        % (1)using Gause's Law J=sigma*E, we can find the current density in both directions:
        Jx = sigma_0*Ex;
        Jy = sigma_0*Ey;
        % now, since we are trying to find the conductivity "seen" by the
        % electrodes, we will use the formula:conductivity = current/potential
        % so, first we will find the current "seen" by the electrodes.
        % we are only interested in the current in the x direction since the
        % potential difference is in the x diraction.
        % (2) calculate the current "seen" by the electrodes.

        I = 1/(h-1)*sum(Jx(YL_idx:YH_idx,h)); %current_seen_by_left_electrode.
        %(the same as the current seen by the right electrode.)
        % (1/h)*1 is the area (dz=1 and dy=1/h)
        %(3)
        conductivity = abs(I)/ delta_V; %units [S/m]
        if Sigma_1 ~= (3)
            disp("Capacity for sigma_1="+sigma_1+": "+C+" [F/["+char(949)+"]m]");
            disp("Conductivity for sigma_1="+sigma_1+": "+conductivity+"[A/V*m]");
        else
            disp("Capacity= "+C+" [F/["+char(949)+"]m]");
            disp("Conductivity= "+conductivity +"[A/V*m]");
        end
    end
  end