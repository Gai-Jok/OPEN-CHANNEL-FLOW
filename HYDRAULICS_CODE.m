%GAI JOK GAI CAID: 44717938
printf("This program was built by GAI JOK[CAID:44717938] \nDISCLAIMER: For any quantity that does not displace units, the program is displaying or requiring the quantity in standard unit.\n")
printf("Flow type:\n1. GRADUALLY VARYING FLOW.\n2. NORMAL FLOW.\n3. RAPIDLY VARYING FLOW.\n");
FLOW = input("Enter the number corresponding to the type of flow: ");

switch FLOW
  case 1
    printf("Solve for:\n1. WIDTH GIVEN.\n2. WIDTH UNKOWN.\n");
choice = input("Enter the number corresponding to if the width is given or not: ");

switch choice
    case 1
% Prompt for input values
S = input("Enter the channel slope (S): ");
n = input("Enter Manning's n: ");
Q = input("Enter flow rate (Q in m^3/s): ");
initial_depth = input("Enter initial depth in meters: ");
target_depth = input("Enter target depth in meters (or 0 to ignore): ");
b = input("Enter the channel width in meters: ");
steps = input("Enter the number of steps for calculation: ");
q = Q / b;

% Constants
g = 9.81; % gravitational acceleration (m/s^2)

% Calculate normal depth (hn) and critical depth (hc) for the flow
hn = ((n * q) / sqrt(S))^(3/5); % Normal depth
hc = ((q^2) / g)^(1/3); % Critical depth

% Determine flow regime
if hn > hc
    flow_type = "subcritical (mild slope)";
    fprintf("The flow is subcritical with a mild slope.\n");
elseif hn < hc
    flow_type = "supercritical (steep slope)";
    fprintf("The flow is supercritical with a steep slope.\n");
else
    flow_type = "critical";
    fprintf("The flow is at critical depth.\n");
end


if initial_depth<target_depth
% Initial values for the iterative GVF calculation
h_delta = (target_depth - initial_depth) / steps; % Step size for distance increment
x = 0; % Initialize total distance along the channel
h = initial_depth; % Start from initial depth

% Iterative GVF calculation
fprintf("Starting GVF calculations...\n");
for i = 1:steps
        h_mid = h + h_delta / 2; % Midpoint depth for this step

        % Calculate slope term m based on midpoint depth
        m = (1 - (q^2 / (g * h_mid^3))) / (S - ((n * q / (h_mid^(5/3)))^2)*((1+((2*h_mid)/b))^(4/3)));

        % Distance increment for this step, ensuring it is positive
        x_delta = abs(m * h_delta);

        % Accumulate total distance
        x = x + x_delta;

        % Update depth for next iteration
        h = h + h_delta;


    % Display current step values
    fprintf("Step %d - Depth: %.4f m, Distance: %.2f m\n", i, h, x);

    % Stop if target depth is reached or exceeded
    if target_depth > 0 && h >= target_depth
        break;
    end
end

% Final result display
fprintf("\nTotal distance from the initial point to reach a depth of %.2f m is approximately %.2f meters using %d steps.\n", h, x, steps);
else
% Set up initial values for GVF calculation
h_delta = (target_depth - initial_depth) / steps; % Step size for distance increment
x = 0; % Initialize total distance along the channel
h = initial_depth; % Start from initial depth
direction = sign(target_depth - initial_depth); % Determine if depth increases (+1) or decreases (-1)

% GVF calculation loop
fprintf("Starting GVF calculations...\n");
for i = 1:steps
        h_mid = h + h_delta / 2; % Midpoint depth for this step
        hs = h + h_delta; % New depth at each step increment

        % Calculate slope term m based on midpoint depth
        m = (1 - (q^2 / (g * h_mid^3))) / (S - ((n * q / (h_mid^(5/3)))^2)*((1+((2*h_mid)/b))^(4/3)));

        % Distance increment for this step, ensuring it is positive
        x_delta = abs(m * h_delta);

        % Accumulate total distance
        x = x + x_delta;

        % Update depth for next iteration
        h = h + h_delta;

    % Display current step values
    fprintf("Step %d - Depth: %.4f m, Distance: %.2f m\n", i, h, x);

    % Stop if target depth is reached or exceeded in the specified direction
    if direction > 0 && h >= target_depth
        break;
    elseif direction < 0 && h <= target_depth
        break;
    end
end

% Final result display
fprintf("\nTotal distance from the initial point to reach a depth of %.2f m is approximately %.2f meters using %d steps.\n", h, x, steps);
endif
case 2
% Prompt for input values
S = input("Enter the channel slope (S): ");
n = input("Enter Manning's n: ");
q = input("Enter flow rate per unit width (q in m^3/s per m): ");
initial_depth = input("Enter initial depth in meters: ");
target_depth = input("Enter target depth in meters (or 0 to ignore): ");
steps = input("Enter the number of steps for calculation: ");

% Constants
g = 9.81; % gravitational acceleration (m/s^2)

% Calculate normal depth (hn) and critical depth (hc) for the flow
hn = ((n * q) / sqrt(S))^(3/5); % Normal depth
hc = ((q^2) / g)^(1/3); % Critical depth

% Determine flow regime
if hn > hc
    flow_type = "subcritical (mild slope)";
    fprintf("The flow is subcritical with a mild slope.\n");
elseif hn < hc
    flow_type = "supercritical (steep slope)";
    fprintf("The flow is supercritical with a steep slope.\n");
else
    flow_type = "critical";
    fprintf("The flow is at critical depth.\n");
end

if initial_depth<target_depth
% Initial values for the iterative GVF calculation
h_delta = (target_depth - initial_depth) / steps; % Depth increment per step
x = 0; % Initialize total distance along the channel
h = initial_depth; % Start from initial depth

% Iterative GVF calculation
fprintf("Starting GVF calculations...\n");
for i = 1:steps
        h_mid = h + h_delta / 2; % Midpoint depth for this step
        hs = h + h_delta; % New depth at each step increment

        % Calculate slope term m based on midpoint depth
        m = (1 - (q^2 / (g * h_mid^3))) / (S - ((n * q / (h_mid^(5/3)))^2));


         % Distance increment for this step, ensuring it is positive
        x_delta = abs(m * h_delta);
       % Accumulate total distance
        x = x + x_delta;

        % Update depth for next iteration
        h = h + h_delta;

        % Display current step values
        fprintf("Step %d - Depth: %.4f m, Distance: %.2f m\n", i, h, x);

        % Stop if target depth is reached or exceeded
        if h >= target_depth

    % Stop if target depth is reached or exceeded
    if target_depth > 0 && h >= target_depth
        break;
    end
end

end
% Final result display
fprintf("\nTotal distance from the initial point to reach a depth of %.2f m is approximately %.2f meters using %d steps.\n", h, x, steps);
else
% Set up initial values for GVF calculation
h_delta = (target_depth - initial_depth) / steps; % Depth increment per step
x = 0; % Initialize total distance along the channel
h = initial_depth; % Start from initial depth
direction = sign(target_depth - initial_depth); % Determine if depth increases (+1) or decreases (-1)

% GVF calculation loop
fprintf("Starting GVF calculations...\n");
for i = 1:steps
        h_mid = h + h_delta / 2; % Midpoint depth for this step
        hs = h + h_delta; % New depth at each step increment

        % Calculate slope term m based on midpoint depth
        m = (1 - (q^2 / (g * h_mid^3))) / (S - ((n * q / (h_mid^(5/3)))^2));


         % Distance increment for this step, ensuring it is positive
        x_delta = abs(m * h_delta);
       % Accumulate total distance
        x = x + x_delta;

        % Update depth for next iteration
        h = h + h_delta;

    % Display current step values
    fprintf("Step %d - Depth: %.4f m, Distance: %.2f m\n", i, h, x);

    % Stop if target depth is reached or exceeded in the specified direction
    if direction > 0 && h >= target_depth
        break;
    elseif direction < 0 && h <= target_depth
        break;
    end
end

% Final result display
fprintf("\nTotal distance from the initial point to reach a depth of %.2f m is approximately %.2f meters using %d steps.\n", h, x, steps);
endif
endswitch
case 2
printf("WHAT IS THE SHAPE OF THE CHANNEL:\n1. RECTANGLE.\n2. TRIANGLE.\n3. TRAPEZIUM.\n4. CIRCLE.\n5. COMPOSITE SHAPE.\n");
shape = input("Enter the number corresponding to the shape of the channel: ");
switch shape
  case 1
    Q=input("Enter the discharge:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
b=input("Enter the channel width:");
g=9.81;
func=@(h) h-((n*Q)/(b*sqrt(S)))^(3/5)*(1+(2*h/b))^(2/5);
h=fzero(func,1);
A=h*b;
V=Q/A;
Fr=V/sqrt(g*h);
q=Q/b;
hc=(q^2/g)^(1/3);
Ec=1.5*hc;
Sc=((n*Q)/b)^(2)*((1+((2*hc)/b))^(4/3))/((hc)^(10/3));
printf("The normal depth is %d.\n",h);
printf("Froude's number is %d.\n",Fr);
printf("The critical depth is %d.\n",hc);
printf("The critical slope is %d.\n",Sc);
printf("The specific energy is %d.\n",Ec);
if h>hc
  printf("The normal depth is subcritical.\n");
else
   printf("The normal depth is subcritical.\n");
   endif
case 2
   Q=input("Enter the discharge:");
x=input("Enter angle on which the tip of triangular is standing:");
h=input("Enter the depth on one side incase of a hydraulic jump: ");
g=9.81;
z=x/2;
hc=((2*(Q^2))/(g*(tand(z))^(2)))^(1/5);
printf("The critical depth is %d.\n",hc);
H=((1/3)*g*(h^3)*tand(z))+((Q^2)/((h^2)*tand(z)))
eq=@(hu) hu-cbrt((3/(g*tand(z)))*(H-((Q^2)/((hu^2)*tand(z)))));
my_intialguess=(3/(g*tand(z)))*H;
hu=fzero(eq,my_intialguess);
printf("The depth on the other side of the jump is %d.\n",hu)
case 3
GIVEN=input("WAS MANNING'S NUMBER(M) OR CHEZY'S CONSTANT(C) GIVEN:","s");
switch GIVEN
  case {"MANNING'S NUMBER","M"}
Q=input("Enter the discharge:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
t=input("Enter channel top width:");
b=input("Enter channel bottom width:");
d=input("Enter channel depth:");
g=9.81;
m=((t-b)/2)/d;
my_func=@(h) ((1/n)*(S^0.5)*((h*(b+(m*h)))/(b+(2*h)*(1+m^2)^0.5))^(2/3)*(h*(b+(m*h))))-Q;
h= fzero(my_func,0.5);
A=h*(b+(m*h));
P=b+((2*h)*sqrt(1+(m)^2));
R=A/P;
V=Q/A;
bs=b+2*m*h;
hm=A/(bs);
Fr=V/sqrt(g*hm);
eqtn=@(hc) (((Q^2)*(b+(2*m*hc)))/(g*(hc*(b+(m*hc)))^3))-1;
hc=fsolve(eqtn,1);
q=Q/b;
func=@(Sc) ((1/n)*(Sc^0.5)*((hc*(b+(m*hc)))/(b+(2*hc)*(1+m^2)^0.5))^(2/3)*(hc*(b+(m*hc))))-Q;
Sc=fzero(func,1);
printf("The normal depth is %d.\n",h);
printf("Froude's number is %d.\n",Fr);
printf("The critical depth is %d.\n",hc);
printf("The critical slope is %d.\n",Sc);
case {"CHEZY'S CONSTANT","C"}
Q=input("Enter the discharge:");
C=input("Enter CHEZY'S CONSTANT:");
S=input("Enter the streamline slope:");
b=input("Enter channel bottom width:");
u=input("Enter the exterior angle of the triangle:");
g=9.81;
m=1/tand(u)
eq=@(h) C*sqrt((h*S*(b+(m*h)))/((2*h*sqrt(1+(m^2)))+b))*h*(b+(m*h))-Q;
h=fsolve(eq,1);
my_eq=@(hc) (((Q^2)*((2*m*hc)+b))/(g*((hc*(b+(m*hc)))^3)))-1;
hc=fsolve(my_eq,1);
printf("The normal depth is %d.\n",h);
printf("The critical depth is %d.\n",hc);
endswitch
case 4
Q=input("Enter the discharge:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
R=input("Enter the radius of the circle:");
g=9.81;
func=@(X) ((1/n)*((((R^2)*(X-(sin(X))*(cos(X))))/(2*R*X))^(2/3))*(S^(0.5))*((R^2)*(X-(sin(X))*(cos(X)))))-Q;
X=fsolve(func,0.5);
A=(R^2)*(X-(sin(X))*(cos(X)));
P=2*R*X;
h=R-(R*cos(X));
bs=2*R*sin(x);
hm=A/bs;
V=Q/A;
Fr=((V)/sqrt(g*hm));
my_func=@(Xc) (((Q^2)*(2*R*sin(Xc)))/(g*((R^2)*(Xc-(sin(Xc))*(cos(Xc))))^(3)))-1;
Xc=fsolve(my_func,1);
hc=R-(R*cos(Xc));
printf("The normal depth is %d.\n",h);
printf("Froude's number is %d.\n",Fr);
printf("The critical depth is %d.\n",hc);
case 5
SHAPE=input("Enter the shape below the section with strainght sides in upper case:","s");
g=9.81;
switch SHAPE
case {"TRIANGLE"}
Q=input("Enter the discharge:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
y=input("Enter the height of triangle:");
x=input("Enter angle on which the tip of triangular is standing:");
g=9.81;
z=x/2;
A=(y^2)*tand(z);
P=(2*y)/cosd(z);
Q_check=((1/n)*((((y^2)*(tand(z)))/((2*y)/cosd(z)))^(2/3))*(S^0.5)*((y^2)*tand(z)));
if Q_check<Q
  eqtn=@(h) ((1/n)*(((y*tand(z)*((2*h)-y))/(((2*(h-y)))+((2*y)/cosd(z))))^(2/3))*(S^0.5)*((y)*(tand(z))*((2*h)-y)))-Q;
  h=fsolve(eqtn,1);
  A=y*tand(z)*((2*h)-y);
  bs=2*y*tand(z);
  hm=A/bs;
  V=Q/A;
  Fr=V/sqrt(g*hm);
  printf("The normal depth is %d.\n",h);
  printf("Froude's number at normal depth is %d.\n",Fr);
else
  my_eqtn=@(h) ((1/n)*((((h^2)*tand(z))/((2*h)/cosd(z)))^(2/3))*(S^0.5)*(h^2)*tand(z))-Q;
  h=fsolve(my_eqtn,0.5);
  A=(h^2)*tand(z);
  bs=2*h*tand(z);
  hm=A/bs;
  V=Q/A;
  Fr=V/sqrt(g*hm);
  printf("The normal depth is %d.\n",h);
  printf("Froude's number at normal depth is %d.\n",Fr);
endif
 func=@(hc) (((Q^2)*(2*hc)*(tand(z)))/((g*((hc^2)*tand(z))^3)))-1;
  hc=fsolve(func,1);
  printf("The critical depth is %d.\n",hc);
case {"SEMICIRCLE"}
SCENARIO=input("Was Manning's number given?(YES or NO):","s");
switch SCENARIO
  case {"NO"}
    Q=input("Enter the discharge:");
    S=input("Enter the streamline slope:");
    R=input("Enter the radius of the circle:");
    d=input("Enter depth of flow:");
    A=0.5*pi*(R^2);
    P=0.5*2*pi*R;
    Rh=A/P;
    n=(1/Q)*(Rh)^(2/3)*(S^0.5)*A;
    printf("Mannings number is %d.\n",n);
    case {"YES"}
    Q=input("Enter the discharge:");
    n=input("Enter manning's number:");
    S=input("Enter the streamline slope:");
    R=input("Enter the radius of the circle:");
    eq=@(h) ((1/n)*((((0.5*pi*(R^2))+(2*R*(h-R)))/((2*(h-R))+(pi*R)))^(2/3))*(S^0.5)*(((0.5*pi*(R^2))+(2*R*(h-R)))))-Q;
    h=fsolve(eq,1)
    A=(0.5*pi*(R^2))+(2*R*(h-R));
    P=(2*h)+(pi*R);
    Rh=A/P;
    V=Q/A;
    bs=2*R;
    hm=A/bs;
    Fr=V/sqrt(g*hm);
if Fr>1
    printf("The normal flow at this discharge is supercritical (Fr > 1); hence, the channel is hydraulically steep.\n");
    printf("The normal depth is %d.\n",h);
    printf("Froude's number is %d.\n",Fr);
 else
    printf("The normal flow at this discharge is subcritical (Fr < 1); hence, the channel is hydraulically mild.\n");
    printf("The normal depth is %d.\n",h);
    printf("Froude's number is %d.\n",Fr);
  endif
endswitch
endswitch
endswitch
case 3
printf("Flow condition:\n1. WEIR.\n2. CONSTRICTION.\n3. EXPANSION.\n4. SLUICEDGATES.\n5. BRIDGE.\n6. BLOCKS.\n7. SPILLWAY.\n");
CONDITION = input("Enter the number corresponding to the flow condition: ");
switch CONDITION
  case 1
    HYDRAULIC=input("IS THE HYDRAULIC JUMP APPRECIATED?(YES OR NO):","s");
switch HYDRAULIC
  case {"NO"}
CONDITION=input("WAS THE WIDTH OF THE CHANNEL GIVEN?(YES OR NO OR DEPRESSION):","s");
switch CONDITION
  case {"NO"}
q=input("Enter the discharge per unit width:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
hx=input("Enter height of the channel(Assume any value greater than zero incase not given):");
z=input("What is the height of the weir:");
g=9.81;
hn=((q*n)/sqrt(S))^(3/5);
Vn=q/hn;
Ea=hn+(Vn^2/(2*g));
hc=((q^2)/g)^(1/3);
Ec=(3/2)*hc;
zs=hx+((q^2)/(2*g*(hx^2)))-(1.5*hc);
Hc=z+Ec;
if Ea>Hc
  E=Ea-z;
  eq=@(h) h-E+((q^2)/(2*g*(h^2)));
  h=fsolve(eq,1);
  x=hn;
  printf("The depth upstream the weir is %d.\n",x);
  printf("The depth downstream the weir is %d.\n",x);
  printf("The depth over the weir is %d.\n",h);
else
  eq=@(h) ((q)/sqrt(2*g*(Hc-h)))-h;
  h=fsolve(eq,1);
  my_eq=@(hk) hk-Hc+((q^2)/(2*g*(hk^2)));
  hk=fsolve(my_eq,1);
  printf("The depth upstream the weir is %d.\n",hk);
  printf("The depth downstream the weir is %d.\n",h);
  printf("The depth over the weir is %d.\n",hc);
endif
zw=Ea-Ec;
printf("The height of the weir to make the flow just critical is %d.\n",zw);
printf("The height of the weir to make the flow overbank is %d.\n",zs);
printf("The normal depth in the channel is %d.\n",hn);
case {"YES"}
Q=input("Enter the discharge:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
b=input("Entee the channel width:");
hx=input("Enter height of the channel:");
z=input("What is the height of the weir(Assume any value greater than zero incase not given):");
g=9.81;
q=Q/b;
func=@(hn) ((1/n)*(((b*hn)/(b+(2*hn)))^(2/3))*sqrt(S)*b*hn)-Q;
hn=fsolve(func,1);
Vn=q/hn;
Ea=hn+(Vn^2/(2*g));
hc=((q^2)/g)^(1/3);
Ec=(3/2)*hc;
Hc=z+Ec;
zs=hx+((q^2)/(2*g*(hx^2)))-(1.5*hc);
if Ea>Hc
  E=Ea-z;
  eq=@(h) h-E+((q^2)/(2*g*(h^2)));
  h=fsolve(eq,1);
  printf("The depth upstream the weir is %d.\n",hn);
  printf("The depth downstream the weir is %d.\n",hn);
  printf("The depth over the weir is %d.\n",h);
else
  eq=@(h) ((q)/sqrt(2*g*(Hc-h)))-h;
  h=fsolve(eq,1);
  my_eq=@(hk) hk-Hc+((q^2)/(2*g*(hk^2)));
  hk=fsolve(my_eq,1);
  printf("The depth upstream the weir is %d.\n",hk);
  printf("The depth downstream the weir is %d.\n",h);
  printf("The depth over the weir is %d.\n",hc);
endif
zw=Ea-Ec;
printf("The height of the weir to make the flow just critical is %d.\n",zw);
printf("The height of the weir to make the flow overbank is %d.\n",zs);
case {"DEPRESSION"}
SCENARIO=input("IS THERE A WEIR IN THE DEPRESSION?(YES OR NO):","s");
switch SCENARIO
  case {"NO"}
q=input("Enter the discharge per unit width:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
z=input("What is the depth of the depression:");
g=9.81;
hn=((q*n)/sqrt(S))^(3/5);
Vn=q/hn;
Ea=hn+(Vn^2/(2*g));
hc=((q^2)/g)^(1/3);
if hn>hc
  my_func=@(ht) Ea+z-((q^2)/(2*g*(ht^2)))-ht;
  ht=fsolve(my_func,1);
  printf("The depth at A is %d.\n",hn)
  printf("The depth at B is %d.\n",hn)
  printf("The depth at C is %d.\n",ht)
else
  disp("SCENARIO NOT CATERED FOR")
endif
case {"YES"}
q=input("Enter the discharge per unit width:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
zu=input("What is the height of the weir:");
z=input("What is the depth of the depression:");
g=9.81;
hn=((q*n)/sqrt(S))^(3/5);
Vn=q/hn;
Ea=hn+(Vn^2/(2*g));
hc=((q^2)/g)^(1/3);
Hc=zu-z+(1.5*hc);
if Hc>Ea
  my_equtn=@(hu) Hc-hu-((q^2)/(2*g*(hu^2)));
  hu=fsolve(my_equtn,1);
  equtn=@(hv) Hc+z-hv-((q^2)/(2*g*(hv^2)));
  hv=fsolve(equtn,1);
  equatn=@(hw) hw-sqrt(((q^2)/(2*g))/(Hc+z-hw));
  hw=fsolve(equatn,1);
  my_equatn=@(ho) ho-sqrt(((q^2)/(2*g))/(Hc-ho));
  ho=fsolve(my_equatn,1);
  printf("The depth at A is %d.\n",hu);
  printf("The depth at B is %d.\n",hv);
  printf("The depth at C is %d.\n",hc);
  printf("The depth at D is %d.\n",hw);
  printf("The depth at E is %d.\n",ho);
else
  disp("SCENARIO NOT CATERED FOR.")
endif
endswitch
endswitch
case {"YES"}
q=input("Enter the discharge per unit width:");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
hx=input("Enter height of the channel(Assume any value greater than zero incase not given):");
z=input("What is the height of the weir:");
g=9.81;
hn=((q*n)/sqrt(S))^(3/5);
Vn=q/hn;
Ea=hn+(Vn^2/(2*g));
hc=((q^2)/g)^(1/3);
Ec=(3/2)*hc;
zs=hx+((q^2)/(2*g*(hx^2)))-(1.5*hc);
Hc=z+Ec;
Fr=Vn/sqrt(g*hn);
zw=Ea-Ec;
my_eq1=@(hm) ((q)/sqrt(2*g*(Hc-hm)))-hm;
hm=fzero(my_eq1,1);
hg=(hn/2)*(sqrt(1+(8*(Fr^2)))-1);
printf("The depth upstream the weir is %d.\n",hg);
printf("The depth downstream the weir is %d.\n",hm);
printf("The height of the weir to make the flow just critical is %d.\n",zw);
printf("The height of the weir to make the flow overbank is %d.\n",zs);
printf("The normal depth in the channel is %d.\n",hn);
printf("The depth over the weir is %d.\n",hc);
if hm>hg
  printf("Since there is no control,the hydraulic jump must actually occur at (or just before) the downstream end of the weir.\n");
else
  printf("Since the is no control,the hydraulic jump would actually occur the downstream end of the weir.\n");
  endif
endswitch
case 2
condition=input("Is the depth given( YES OR NO): ", "s");
switch condition
  case {"NO"}
Q = input("Enter the discharge (Q) in cubic meters per second: ");
n = input("Enter Manning's roughness coefficient (n): ");
S = input("Enter the channel slope (S): ");
b = input("Enter the initial channel width (b): ");
bm = input("Enter the channel width at constriction (bm): ");
z = input("Enter the elevation change (+ for rise, - for drop): ");
g = 9.81;

func = @(ha) ha - ((n * Q) / (b * sqrt(S)))^(3/5) * (1 + (2 * ha / b))^(2/5);
ha = fzero(func, 6);
A = ha * b;
V = Q / A;
Fr = V / sqrt(g * ha);
q = Q / b;
hc = (q^2 / g)^(1/3);
Ec = 1.5 * hc;
Sc = ((n * Q) / b)^2 * (1 + (2 * hc / b))^(4/3) / hc^(10/3);
Ha = ha + (V^2) / (2 * g)
Hc = z + Ec;
qm = Q / bm;
hcm = (qm^2 / g)^(1/3);
Ecm = 1.5 * hcm;
Hcm = z + Ecm
A_constriction = hcm * bm;
V_constriction = Q / A_constriction;
Fr_constriction = V_constriction / sqrt(g * htm)

is_critical = (abs(Fr_constriction - 1) < 1e-3)

if (z == 0) && (is_critical == 0)
    my_eq = @(h) Hcm - h - (Q^2 / (2 * g * bm^2 * h^2));
    h = fsolve(my_eq, 5);
    printf("The depth upstream is %.6f.\n", ha);
    printf("The depth downstream is %.6f.\n", ha);
    printf("The depth at the throat is %.6f.\n", h);
elseif (z == 0) && (is_critical == 1)
    my_eq = @(h) Hcm - h - (Q^2 / (2 * g * b^2 * h^2));
    h = fsolve(my_eq, 5);
    eq = @(hu) hu - (sqrt(((Q^2) / (2 * g * b^2)) /(Hcm - hu)));
    hu = fzero(eq, 0);
    printf("The depth upstream is %.6f.\n", h);
    printf("The depth downstream is %.6f.\n", hu);
    printf("The depth at the throat is %.6f.\n", hc);
    Vx= Q / ( b * h );
    Vu= Q / ( b * hu );
    Vt= Q / ( b * hc );
    printf("The velocity upstream is %.6f.\n", Vx);
    printf("The velocity downstream is %.6f.\n", Vu);
    printf("The velocity at the throat is %.6f.\n", Vt);
elseif (z > 0)
    my_eq = @(hr) Hcm - hr - (Q^2 / (2 * g * b^2 * hr^2));
    hr = fsolve(my_eq, 4);
    eq = @(hu) hu - (sqrt(((Q^2) / (2 * g * b^2)) /(Hcm - hu)));
    hu = fzero(eq, 0);
    printf("The depth upstream is %.6f.\n", hr);
    printf("The depth downstream is %.6f.\n", hu);
    printf("The depth at the throat is %.6f.\n", hcm);

else  % This is the case where z < 0
    my_eq = @(hp) Ha - z - hp - (Q^2 / (2 * g * bm^2 * hp^2));
    hp = fsolve(my_eq, 6);
    printf("The depth upstream is %.6f.\n", ha);
    printf("The depth downstream is %.6f.\n", ha);
    printf("The depth at the throat is %.6f.\n", hp);
endif

printf("The normal depth is %.6f.\n", ha);
printf("Froude's number is %.6f.\n", Fr);
printf("The critical depth is %.6f.\n", hc);
printf("The critical slope is %.6f.\n", Sc);
printf("The specific energy is %.6f.\n", Ec);
printf("The critical depth at the constriction is %.6f.\n", hcm);
printf("The specific energy at the constriction is %.6f.\n", Ecm);

if is_critical
    printf("Flow becomes critical at the constriction.\n");
else
    printf("Flow does not reach critical conditions at the constriction.\n");
endif

zb = Ha-Ecm;
     printf("The bed must be raised to %d so that the total head under critical conditions equals that in normal flow.\n",zb);
   case {"YES"}
Q = input("Enter the discharge (Q) in cubic meters per second: ");
n = input("Enter Manning's roughness coefficient (n): ");
S = input("Enter the channel slope (S): ");
b = input("Enter the initial channel width (b): ");
bm = input("Enter the channel width at constriction (bm): ");
z = input("Enter the elevation change (+ for rise, - for drop): ");
hacm= input("Enter the depth of the channel: ");
g = 9.81;

A = hacm * bm
V = Q / A
Fr = V / sqrt(g * hacm)
q = Q / b
hc = (q^2 / g)^(1/3)
Ec = 1.5 * hc
Sc = ((n * Q) / b)^2 * (1 + (2 * hc / b))^(4/3) / hc^(10/3)
Hacm = hacm + (V^2) / (2 * g)
Hc = z + Ec
qm = Q / bm
hcm = (qm^2 / g)^(1/3)
Ecm = 1.5 * hcm;
Hcm = z + Ecm
A_constriction = hacm * bm;
V_constriction = Q / A_constriction;
Fr_constriction = V_constriction / sqrt(g * hacm)

is_critical = (abs(Fr_constriction - 1) < 1e-3)

if (z == 0) && (is_critical == 0)
    my_eq = @(h) Hacm - h - (Q^2 / (2 * g * b^2 * h^2));
    h = fsolve(my_eq, 5)
    printf("The depth upstream is %.6f.\n", h);
    printf("The depth downstream is %.6f.\n", h);
    Vx= Q / ( b * h )
    Vu= Q / ( b * h )
    printf("The velocity upstream is %.6f.\n", Vx);
    printf("The velocity downstream is %.6f.\n", Vu);
elseif (z == 0) && (is_critical == 1)
    my_eq = @(h) Hacm - h - (Q^2 / (2 * g * b^2 * h^2));
    h = fsolve(my_eq, 5);
    eq = @(hu) hu - (sqrt(((Q^2) / (2 * g * b^2)) /(Hacm - hu)));
    hu = fzero(eq, 0);
    printf("The depth upstream is %.6f.\n", h);
    printf("The depth downstream is %.6f.\n", hu);
    Vx= Q / ( b * h );
    Vu= Q / ( b * hu );
    printf("The velocity upstream is %.6f.\n", Vx);
    printf("The velocity downstream is %.6f.\n", Vu);
elseif (z > 0)
    my_eq = @(hr) Hacm - hr - (Q^2 / (2 * g * b^2 * hr^2));
    hr = fsolve(my_eq, 4);
    eq = @(hu) hu - (sqrt(((Q^2) / (2 * g * b^2)) /(Hacm - hu)));
    hu = fzero(eq, 0);
    printf("The depth upstream is %.6f.\n", hr);#
    printf("The depth downstream is %.6f.\n", hu);
    printf("The depth at the throat is %.6f.\n", hcm);

else  % This is the case where z < 0
    my_eq = @(hp) Hacm - z - hp - (Q^2 / (2 * g * bm^2 * hp^2));
    hp = fsolve(my_eq, 6);
    printf("The depth upstream is %.6f.\n", ha);
    printf("The depth downstream is %.6f.\n", ha);
    printf("The depth at the throat is %.6f.\n", hp);
endif

printf("Froude's number is %.6f.\n", Fr);
endswitch
case 3
disp("Assumptions:")
disp("● The hydraulic jump is triggered immediately at the expansion.")
disp("● Reactions from the expansion end walls are in equilibrium with a hydrostatic pressure distribution.")
Q = input("Enter the discharge (Q) in cubic meters per second: ");
b = input("Enter the initial channel width (b): ");
B = input("Enter the channel width after the expansion (B): ");
hacm= input("Enter the depth upstream of the channel: ");
g = 9.81;
rho = 1000;

V = Q / (b * hacm);
my_func=@(hw) hw-sqrt((2/(rho*g*B))*((rho*Q*V)+(0.5*rho*g*(hacm^2)*B)-(((rho*(Q^2))/(B*hw)))));
my_intialguess=((2/(rho*g*B))*((rho*Q*V)+(0.5*rho*g*(hacm^2)*B)));
hw=fzero(my_func,my_intialguess);
printf("The depth downstream is %dm.\n",hw)
case 4
printf("Solve for:\n1. Only downstream depth given.\n2. Both downstream and upstream depths given.\n3.Only upstream depth given.\n");
choice = input("Enter the number corresponding to the known variable: ");

switch choice
    case 1
Q = input("Enter the discharge (Q) in cubic meters per second: ");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
b = input("Enter the channel width (b): ");
h = input("Enter the depth downstream: ");
g = 9.81;
rho = 1000;

 V_downstream  = Q / (b * h);
 H_t = h + ( V_downstream  ^2)/(2*g);
 my_eq = @(hr) H_t - hr - (Q^2 / (2 * g * b^2 * hr^2));
 hr = fsolve(my_eq, 4);
 V_upstream = Q / (b * hr);
 Fr_upstream= V_upstream / sqrt( g * hr);
 Fr_downstream= V_downstream / sqrt( g * h);
 F=(0.5*rho*g*b*((hr^2)-(h^2)))-(rho*((Q^2)/b)*((1/h)-(1/hr)));
 func=@(ha) ha-((n*Q)/(b*sqrt(S)))^(3/5)*(1+(2*ha/b))^(2/5);
ha=fzero(func,1);
q=Q/b;
hc=(q^2/g)^(1/3);
Ec=1.5*hc;
Sc=((n*Q)/b)^(2)*((1+((2*hc)/b))^(4/3))/((hc)^(10/3));
printf("The normal depth is %d.\n",ha);
printf("The critical depth is %d.\n",hc);
printf("The critical slope is %d.\n",Sc);
 printf("The total head given downstream depth is %dm.\n",H_t)
 printf("The depth upstream is %dm.\n",hr)
 printf("The velocity of flow upstream is %d m/s.\n",V_upstream)
 printf("Froude's number upstream is %d.\n",Fr_upstream)
 printf("Froude's number downstream is %d.\n",Fr_downstream)
 printf("The force on the gate is %d.\n",F)
 if (ha<hc)
   if hr>hc
   printf("The hydrualic jump will occur upstream.\n")
   q = Q / b;
Vu = q / ha;
 Fr_downstream= Vu / sqrt( g * ha);
 hk=(ha/2)*(sqrt(1+(8*(Fr_downstream^2)))-1);
 printf("The depth downstream the hydraulic jump is %d.\n",hk)
 printf("The depth upstream the hydraulic jump is %d.\n",ha)
 elseif hr<hc
    printf("The hydrualic jump will not occur upstream.\n")
    elseif h>hc
   printf("The hydrualic jump will occur downstream.\n")
 elseif h<hc
    printf("The hydrualic jump will not occur downstream.\n")
    endif
 else
    if hr>hc
   printf("The hydrualic jump will not occur upstream.\n")
 elseif hr<hc
    printf("The hydrualic jump will occur upstream.\n")
    elseif h>hc
   printf("The hydrualic jump will not occur downstream.\n")
 elseif h<hc
    printf("The hydrualic jump will occur downstream.\n")
     endif
     endif
case 2
 b = input("Enter the channel width (b): ");
h = input("Enter the depth downstream: ");
hu = input("Enter the depth upstream: ");
n=input("Enter manning's number:");
S=input("Enter the streamline slope:");
g = 9.81;
rho = 1000;

my_eq=@(q) hu-h-(((q^2)/(2*g))*((1/(h^2))-(1/(hu^2))));
q=fsolve(my_eq,0);
Q=b*q;
func=@(hn) hn-((n*Q)/(b*sqrt(S)))^(3/5)*(1+(2*hn/b))^(2/5);
hn=fzero(func,1);
Vu = q / hu;
V = q / h;
 F=(0.5*rho*g*b*((hu^2)-(h^2)))-(rho*q*b*((V-Vu)));
 Fr_downstream= V / sqrt( g * h);
 hk=(h/2)*(sqrt(1+(8*(Fr_downstream^2)))-1);
 Vk = q / hk;
 Ha=h+((V^2)/(2*g));
 Hb=hk+((Vk^2)/(2*g));
 x=((Ha-Hb)/Ha);
 Fn=(0.5*rho*g*b*((h^2)-(hn^2)))+(rho*q*b*(V-(q/hn)));
 printf("The discharge is %d.\n",Q)
 printf("The normal depth at this discharge is %d.\n",hn);
  printf("The force on the gate is %d.\n",F)
   printf("The depth downstream the hydraulic jump is %d.\n",hk)
    printf("The fraction of the fluid energy that is dissipated in the jump is %d.\n",x)
      printf("If a block is placed on the downstream to cause a hydraulic jump, the force on the blocks is %d.\n",Fn)
  case 3
Q = input("Enter the discharge (Q) in cubic meters per second: ");
b = input("Enter the channel width (b): ");
hu = input("Enter the depth upstream: ");
g = 9.81;
rho = 1000;

condition=input("Does an energy loss occur.(YES or NO):  ","s");
switch condition
  case {"NO"}
 V_upstream  = Q / (b * hu);
 H_tm = hu + ( V_upstream  ^2)/(2*g);
    eq = @(hs) hs - (sqrt(((Q^2) / (2 * g * b^2)) /(H_tm - hs)));
    hs = fzero(eq, 0);
    V_downstream = Q / ( b * hs);
 Fr_downstream=  V_downstream / sqrt( g * hs);
 printf("The depth downstream is %d.\n",hs)
 printf("Froude's number downstream is %d.\n",Fr_downstream)
case {"YES"}
 k=input("Enter the percentage of energy lost:  ");
 V_upstream  = Q / (b * hu);
 H_tm = hu + ( V_upstream  ^2)/(2*g);
 H_tmk=((100-k)/100)*H_tm
    eq = @(hs) hs - (sqrt(((Q^2) / (2 * g * b^2)) /(H_tmk - hs)));
    hs = fzero(eq, 0);
    V_downstream = Q / ( b * hs);
 Fr_downstream=  V_downstream / sqrt( g * hs);
 printf("The depth downstream is %d.\n",hs)
 printf("Froude's number downstream is %d.\n",Fr_downstream)
  endswitch
endswitch
case 5
HYDRAULIC=input("IS THE HYDRAULIC JUMP APPRECIATED?(YES OR NO):","s");
switch HYDRAULIC
  case {"NO"}
disp("Assumptions:")
disp("● constant total head (main assumption).")
disp("● constricted section long enough to establish parallel flow with critical depth.")
disp("● downstream controls do not prevent supercritical flow being established.")
Q = input("Enter the discharge (Q) in cubic meters per second: ");
b = input("Enter the initial channel width (b): ");
bm = input("Enter the channel width at constriction (bm): ");
z=input("Enter height between the bottom of the bridge and the bed of the river: ");
g = 9.81;

qm = Q / bm;
hcm = (qm^2 / g)^(1/3);
Ecm = 1.5 * hcm;
    my_eq = @(hr) Ecm - hr - (Q^2 / (2 * g * b^2 * hr^2));
    hr = fsolve(my_eq, 4);
    eq = @(hu) hu - (sqrt(((Q^2) / (2 * g * b^2)) /(Ecm - hu)));
    hu = fzero(eq, 0);
    printf("The depth upstream the bridge is %.6f.\n", hr);
    printf("The depth downstream the bridge is %.6f.\n", hu);
    printf("The depth under the bridge is %.6f.\n", hcm);
      if hcm>z
      printf("This exceeds the clearance of the bridge deck(%d). Hence, critical conditions cannot be attained and the flow must be choked.\n", z)
    else
      printf("This is less than the clearance of the bridge deck (%d). Hence, critical conditions can be attained and the flow will not be choked.\n", z)
    endif
  case {"YES"}
    disp("Assumptions:")
disp("● constant total head (main assumption).")
disp("● constricted section long enough to establish parallel flow with critical depth.")
disp("● downstream controls do not prevent supercritical flow being established.")
Q = input("Enter the discharge (Q) in cubic meters per second: ");
b = input("Enter the initial channel width (b): ");
bm = input("Enter the channel width at constriction (bm): ");
z=input("Enter height between the bottom of the bridge and the bed of the river: ");
g = 9.81;

qm = Q / bm;
hcm = (qm^2 / g)^(1/3);
Ecm = 1.5 * hcm;
    my_eq = @(hr) Ecm - hr - (Q^2 / (2 * g * b^2 * hr^2));
    hr = fsolve(my_eq, 4);
    eq = @(hu) hu - (sqrt(((Q^2) / (2 * g * b^2)) /(Ecm - hu)));
    hu = fzero(eq, 0);
    V_downstream= Q / ( b * hu);
    Fr_downstream= V_downstream / sqrt( g * hu);
    hg=(hu/2)*(sqrt(1+(8*(Fr_downstream^2)))-1);
    V_upstream= Q / ( b * hu);
    Fr_upstream= V_upstream / sqrt( g * hu);
    hk=(hr/2)*(sqrt(1+(8*(Fr_upstream^2)))-1);
    printf("The depth upstream is %d.\n",hg);
    printf("The depth downstream is %d.\n",hk);
    if hcm>z
      printf("This exceeds the clearance of the bridge deck(%d). Hence, critical conditions cannot be attained and the flow must be choked.\n", z)
    else
      printf("This is less than the clearance of the bridge deck (%d). Hence, critical conditions can be attained and the flow will not be choked.\n", z)
    endif
endswitch
case 6
Q = input ("Enter the discharge: ");
b = input ("Enter the width of the channel: ");
hb = input ("Enter the height of the blocks: ");
hk = input ("Enter the depth upstream: ");
Cd = input ("Enter the drag coefficient: ");
r = input ("Enter the number of blocks: ");
g = 9.81;
rho = 1000;

V = Q / ( b * hk );
F = r * Cd * 0.5 * rho * (V^2) * hb * b;
func=@(hu) (rho*Q*V*(hk/hu))+(0.5*rho*g*(hu^2)*b)+F-(rho*Q*V)-(0.5*rho*g*(hk^2)*b);
hu=fsolve(func,1);
my_func=@(hw) hw-sqrt((2/(rho*g*b))*((rho*Q*V)+(0.5*rho*g*(hk^2)*b)-F-((rho*Q*V*hk)/hw)));
my_intialguess=((2/(rho*g*b))*((rho*Q*V)+(0.5*rho*g*(hk^2)*b)))-F;
hw=fzero(my_func,my_intialguess);
printf("If the hydraulic jump occurs, the depth downstream is %d.\n",hw);
printf("If the hydraulic jump does not occur, the depth downstream is %d.\n",hu);
case 7
   H=input("Enter the total head:");
   n=input("Enter manning's number:");
   S=input("Enter the streamline slope:");
   Sx=input("Enter the new streamline slope(brought about by adding foreign material in the channel):");
   b=input("Enter the channel width:");
   g=9.81;
   hc=(2/3)*H;
   q=(g^0.5)*(hc^1.5);
   Q=q*b;
   func=@(h) h-((n*Q)/(b*sqrt(S)))^(3/5)*(1+(2*h/b))^(2/5);
   h=fzero(func,1);
   func=@(hx) hx-((n*Q)/(b*sqrt(Sx)))^(3/5)*(1+(2*hx/b))^(2/5);
   hx=fzero(func,1);
   F=(0.5*rho*g*b*((h^2)-(hx^2)))+(rho*q*b*(((q/h)-(q/hx))));
   disp("Assuming the flow over the spillway is critical,");
   printf("The normal depth at this discharge is %d.\n",h);
   printf("The discharge is %d.\n",Q);
   printf("The force on the gate is %d.\n",F)
   printf("The normal depth downstream after adding the blocks is %d.\n",hx);
if h>hc
  printf("The slope is hydraulically steep at this discharge.\n");
else
   printf("The slope is not hydraulically steep at this discharge.\n");
 endif
endswitch
endswitch
