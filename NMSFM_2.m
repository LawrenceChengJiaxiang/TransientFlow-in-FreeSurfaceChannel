    %% Calculate the initial condition
    clear all;

    global Q1 Q2 dtrec

    b = 2.5; % Canal width
    g = 9.8;
    n = 0.029; % Manning coefficient

    t = 0;
    tf = 1;
    dt = 0.01;
    dtrec = 3600;
    alpha = 0.6;

    L = 7000; % Domain length
    N = 71; % Number of points in the mesh
    I = 0.0003; % Slope
    h = 0.7;
    dx = L/(N-1); % Space step

    t1 = 2880;
    t2 = 25920;
    Q0 = [0.2 0.2 0.2];
    Q1 = 0.2; % Upstream flow
    Q2 = 0.89;

    Qunit = Q1/b; % Unit flow
    Hn = (Qunit^2 * n^2/I)^0.3; % Normal water depth
    Hc = (Qunit^2/g)^(1/3); % Critical water depth

    h0(N) = h;
    for i = N-1:-1:1
        h0(i) = h0(i+1) - dx*I*(1-(Hn/h0(i+1))^(10/3))/(1-(Hc/h0(i+1))^3)
    end
    for i = 1:N
        u0(i) = Q1/(b*h0(i));
    end

    % figure(1)
    % subplot(2,1,1);plot(0:100:7000,h0);
    % subplot(2,1,2);plot(0:100:7000,u0);

    %% Test the LAX Scheme
    % h1 = h0; u1 = u0;
    % while (t < tf)
    %     t = t + dt;
    %     for i = 2 : N-1
    %         J0(i) = n^2 * u0(i)^2/ h0(i)^(4/3);
    %         h1(i) = dt*(h0(i)*((u0(i-1)-u0(i+1))/(2*dx))+u0(i)*((h0(i-1)-...
    %             h0(i+1))/(2*dx)))+alpha*h0(i)+(1-alpha)*((h0(i+1)+h0(i-1))/2);
    %         u1(i) = g*dt*(I-J0(i)+u0(i)*((u0(i-1)-u0(i+1))/(2*g*dx))+(h0(i-1)-...
    %             h0(i+1))/(2*dx))+alpha*u0(i)+(1-alpha)*((u0(i+1)+u0(i-1))/2);
    %     end
    % end
    % 
    % figure(2)
    % subplot(2,1,1);plot(0:100:7000,h1);
    % subplot(2,1,2);plot(0:100:7000,u1);

    %% Iteration Simulation
    x0 = 0;
    xN = 7000;
    Up = u0; Ur = u0;
    Hp = h0; Hr = h0;
    dk = 0.1/dtrec;
    for i = 1 : N
        Jp(i) = n^2 * Up(i)^2 / Hp(i)^(4/3);
    end

    % figure(3) % This is for the static figure
    for k = 0 : dk : 12

        if mod(k,100/dtrec)==0
            % Release this one if you want the static figure
            % subplot(1,2,1); plot(0:100:7000,Hp(1:71)); hold on
            % subplot(1,2,2); plot(0:100:7000,Up(1:71)); hold on
            kn = uint32(k*dtrec);
            Hmovie(kn+1,1:71) = Hp(1:71);
            Umovie(kn+1,1:71) = Up(1:71);
        end
    %     if ((k >= 0) && (k <= 4))
    %         if mod(k,200/dtrec)==0
    %             subplot(1,2,1);plot(0:100:7000,Hp(1:71))
    %             hold on;
    %             subplot(1,2,2);plot(0:100:7000,Up(1:71))
    %             hold on;
    %         end
    %     else
    %         if mod(k,2000/dtrec)==0
    %             subplot(1,2,1);plot(0:100:7000,Hp(1:71))
    %             hold on;
    %             subplot(1,2,2);plot(0:100:7000,Up(1:71))
    %             hold on;
    %         end
    %     end

        Qi1 = Qinput(k);
        Qi2 = Qinput(k+dk);

        A1 = [1+g*dt*Jp(1)/Up(1), -sqrt(g/Hp(1))-2*g*dt*Jp(1)/(3*Hp(1)),... 
            -1+g*dt*Jp(1)/Up(1), sqrt(g/Hp(1))-2*g*dt*Jp(1)/(3*Hp(1)), 0;
            1/Up(1), 1/Hp(1), 0, 0, 0;
            0, 0, dx, 0, Up(1)-Up(2);
            0, 0, 0, dx, Hp(1)-Hp(2);
            1/2, -sqrt(g/Hp(1))/4, 1/2, -sqrt(g/Hp(1))/4, 1/dt];

        B1 = [g*dt*(I-Jp(1)/3);
            1+Qi2/Qi1;
            Up(1)*dx+x0*(Up(1)-Up(2));
            Hp(1)*dx+x0*(Hp(1)-Hp(2));
            sqrt(g*Hp(1))/2+x0/dt];

        X1 = PG(A1, B1);
        Ur(1) = X1(1); % obtain the current velocity at position 0
        Hr(1) = X1(2);

        An = [1+g*dt*Jp(N)/Up(N), -1+g*dt*Jp(N)/Up(N), -sqrt(g/Hp(N))-2*...
            g*dt*Jp(N)/(3*Hp(N)), 0;
            1/2, 1/2, sqrt(g/Hp(N))/4, 1/dt;
            0, 0, dx, Hp(N-1)-Hp(N);
            0, dx, 0, Up(N-1)-Up(N)];

        Bn = [g*dt*(I-Jp(N)/3)-sqrt(g*Hp(N));
            -(3/4)*sqrt(g*Hp(N))+xN/dt;
            Hp(N)*dx+xN*(Hp(N-1)-Hp(N));
            Up(N)*dx+xN*(Up(N-1)-Up(N))];

        Xn = PG(An, Bn);
        Ur(N) = Xn(1); % obtain the current velocity at position N
        Hr(N) = h; % the height at position N is always h

        for i = 1 : N
            Jr(i) = n^2 * Ur(i)^2 / Hr(i)^(4/3);
        end

        for i = 2 : N-1
            Hr(i) = dt*(Hr(i)*((Ur(i-1)-Ur(i+1))/(2*dx))+Ur(i)*((Hr(i-1)-...
                Hr(i+1))/(2*dx)))+alpha*Hr(i)+(1-alpha)*((Hr(i+1)+Hr(i-1))/2);
            Ur(i) = g*dt*(I-Jr(i)+Ur(i)*((Ur(i-1)-Ur(i+1))/(2*g*dx))+(Hr(i-1)-...
                Hr(i+1))/(2*dx))+alpha*Ur(i)+(1-alpha)*((Ur(i+1)+Ur(i-1))/2);
        end

        Hp = Hr;
        Up = Ur;
        Jp = Jr;

    end

    %% Movie Generator
    figure(4)
    [rowH, colH] = size(Hmovie);
    sceneNum = (rowH-1)/100;
    M = moviein(sceneNum); % Initialize the movie matrix
    Qin = []; % Initialize the sequence of the input flow

    for k = 1:sceneNum % Generate each scene for the movie
       Qin = [Qin, Qinput(k*12/sceneNum)]; 
       subplot(3,1,1);plot(0:12/sceneNum:12*(k-1)/sceneNum,Qin(1:k));
       xlim([0 12]);ylim([0 1]); title('Input volume flow');
       xlabel('time/h');ylabel('volume/(m^3/s)');
       subplot(3,1,2);plot(0:100:7000,Hmovie(k*100+1,1:71));
       xlim([0 7000]);ylim([0.25 1]); title('Height of water downstream');
       xlabel('position length/m');ylabel('height/m');
       subplot(3,1,3);plot(0:100:7000,Umovie(k*100+1,1:71));
       xlim([0 7000]);ylim([0.25 1]); title('Average velocity of water downstream');
       xlabel('position length/m');ylabel('average velocity/(m/s)');
       % Call getframe to transfer each plot to scene
       M(k) = getframe;
    end
    % Call movie function to play the scenes for 1 time.
    movie(M, 1);
