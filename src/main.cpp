#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
//#include "matplotlibcpp.h"
#include <chrono>

//DONE: change applied_latency to 100 for submission
int applied_latency = 100; //latency to apply before sending command, in milliseconds

//double v_target = 10; //speed target
//double T = 3.0; //time ahead to predict
//double dt = 0.1; //timestep to use for prediction

//double v_target = 30;
//double T = 3.0; 
//double dt = 0.1;

// double v_target = 30;//predicted path stretches beyond input path
// double T = 1.0; 
// double dt = 0.05; 
// double v_target = 30; //still too far; reduce again
// double T =  0.5;  //time ahead to predict
// double dt = 0.02; //timestep to use for prediction

//double v_target = 30; //ERROR: was looking at poly track instead of MPC when I said predicted path stretches too far
//double T =  0.2;  //time ahead to predict
//double dt = 0.02; //timestep to use for prediction

// double v_target = 50; //change to ms for int math
// int T   =  1000;      //time ahead to predict, in ms
// int dtm =    50;      //timestep to use for prediction, in ms

// double v_target = 50; //looks good; try higher resolution, with eye on increasing effective_latency
// int T   =  1000;      //time ahead to predict, in ms
// int dtm =    20;      //timestep to use for prediction, in ms

// double v_target = 50; //previous was way too slow; will reduce size significantly
// int T   =   500;      //time ahead to predict, in ms
// int dtm =    50;      //timestep to use for prediction, in ms
// int effective_latency = applied_latency + 150; //account for additional latency of this program: with current settings, measure 100ms in computation alone

// double v_target = 50; //adding effective latency made it worse; try reducing
// int T   =   500;      //time ahead to predict, in ms
// int dtm =    50;      //timestep to use for prediction, in ms
// int effective_latency = applied_latency + 50; 

// double v_target = 50; //my mistake earlier- computation is only 10ms; 100ms is full cycle latency (which should still add to applied latency, but I don't know when in that 100ms the new value actually gets applied
// int T   =   500;      //time ahead to predict, in ms
// int dtm =    50;      //timestep to use for prediction, in ms
// int effective_latency = applied_latency + 0; 

// double v_target = 50; //better, but doesn't look far enough ahead (ideally, would like to match input data); increase T
// int T   =  1000;      //time ahead to predict, in ms
// int dtm =    50;      //timestep to use for prediction, in ms
// int effective_latency = applied_latency + 0; 

double v_target = 50; //MUCH BETTER: will keep this and tune cost function
int T   =  1000;      //time ahead to predict, in ms
int dtm =    50;      //timestep to use for prediction, in ms
int effective_latency = applied_latency + 0; 


size_t N = size_t(T/dtm); //number of timesteps to predict over
double dt = (double)dtm/1000.0;

size_t n_latency = size_t(effective_latency/dtm);



// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  double ecks = 1;
  for (int i = 0; i < coeffs.size(); i++) {
    //result += coeffs[i] * pow(x, i);
    result += coeffs[i] * ecks; //pow is expensive
    ecks *= x;
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

Eigen::VectorXd rel_coord_transform(Eigen::VectorXd abs_coord, Eigen::VectorXd cur_coord, double cur_psi){
  Eigen::VectorXd rel_coord(2);
  Eigen::VectorXd ro_coord(2);
  //double neg_psi = -1 * cur_psi;

  rel_coord[0] = abs_coord[0] - cur_coord[0];
  rel_coord[1] = abs_coord[1] - cur_coord[1];
  
  //ro_coord[0] = rel_coord[0] * cos(neg_psi) - rel_coord[1] * sin(neg_psi);
  //ro_coord[1] = rel_coord[1] * cos(neg_psi) + rel_coord[0] * sin(neg_psi);
  ro_coord[0] = rel_coord[0] * cos(cur_psi) + rel_coord[1] * sin(cur_psi);
  ro_coord[1] = rel_coord[1] * cos(cur_psi) - rel_coord[0] * sin(cur_psi);

  return ro_coord;
}

int main() {
  uWS::Hub h;

  cout << "T: "<<T<<" ms \n";
  // MPC is initialized here!
  MPC mpc;
  mpc.init(N,dt,n_latency,v_target);
  //auto start_time = std::chrono::system_clock::now();
  //h.onMessage([&start_time,&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          // vector<double> ptsx_vec = j[1]["ptsx"];
          // vector<double> ptsy_vec = j[1]["ptsy"];
          // Eigen::VectorXd ptsx(ptsx_vec.size());
          // Eigen::VectorXd ptsy(ptsy_vec.size());
          // for (uint i = 0; i < ptsx_vec.size(); i++){
          //   ptsx[i] = ptsx_vec[i];
          //   ptsy[i] = ptsy_vec[i];
          // }
          vector<double> raw_ptsx = j[1]["ptsx"];
          vector<double> raw_ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          //don't use these for latency- instead use the predictive model I already have!
          //double psteer = j[1]["steering_angle"];
          //double pthrottle = j[1]["throttle"];
          //auto cur_time = std::chrono::system_clock::now();
          //std::chrono::duration<double> elapsed_seconds = cur_time - start_time;
          //cout <<"time "<<elapsed_seconds.count()<<" : car state is x,y,psi,v,psteer,pthrottle "<<px<<","<<py<<","<<psi<<","<<v<<","<<psteer<<","<<pthrottle<<"\n";

          Eigen::VectorXd rel_coord(2);
          Eigen::VectorXd abs_coord(2);
          Eigen::VectorXd cur_coord(2);
          cur_coord[0] = px;
          cur_coord[1] = py;

          //convert reference line into car's coordinates
          vector<double> ref_x_vals;
          vector<double> ref_y_vals;
          Eigen::VectorXd ptsx(raw_ptsx.size());
          Eigen::VectorXd ptsy(raw_ptsy.size());

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          //cout <<"reference points, transformed from point: ("<<cur_coord[0]<<","<<cur_coord[1]<<") and angle "<<psi<<"\n";
          uint pts_index = 0;
          for (uint i = 0; i < raw_ptsx.size(); i++){
            abs_coord[0] = raw_ptsx[i];
            abs_coord[1] = raw_ptsy[i];
            rel_coord = rel_coord_transform(abs_coord,cur_coord,psi);
            //bad idea: same problem, just shifted.  Will control by minimizing steering angle changes instead
            // if(rel_coord[0] < 0){//skip points that are behind car- don't want them to swing plotted path
            //   ptsx.resize(ptsx.size() - 1);
            //   ptsy.resize(ptsy.size() - 1);
            //   continue;
            // }
            ptsx[pts_index] = rel_coord[0];
            ptsy[pts_index] = rel_coord[1];
            ref_x_vals.push_back(rel_coord[0]);
            ref_y_vals.push_back(rel_coord[1]);
            pts_index++;
            //cout << "\t("<<abs_coord[0]<<","<<abs_coord[1]<<")\t("<<rel_coord[0]<<","<<rel_coord[1]<<")\n";
          }
          //cout <<"calling polyfit\n";
          //auto coeffs = polyfit(ptsx,ptsy,3);//create a 3rd order polynomial to plot the track in front
          auto coeffs = polyfit(ptsx,ptsy,3);//now using car's frame of reference, rather than absolute
          //auto coeffs = polyfit(ptsx,ptsy,2);//try lower order poly, since I'm not looking that far ahead
          
          
          //compute cte.  Assume that first point listed in ref_x_vals is where we should be.
          //double delta_x = ptsx[1] - ptsx[0];
          //double delta_y = ptsy[1] - ptsy[0];
          //cte formula from https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
          //double cte = fabs(delta_y*px - delta_x*py + ptsx[1]*ptsy[0] - ptsy[1]*ptsx[0]) / pow(pow(delta_y,2) + pow(delta_x,2),0.5);
          //in cars coordnates, px and py are zero, so terms disappear
          //double cte = fabs(ptsx[1]*ptsy[0] - ptsy[1]*ptsx[0]) / pow(pow(delta_y,2) + pow(delta_x,2),0.5);
          double cte = polyeval(coeffs,0);
          
          //double epsi = psi - atan2(delta_y, fabs(delta_x) > 1e-6 ? delta_x : 1e-6);//should have normalized after subtraction
          //double epsi = 0 - atan2(delta_y, fabs(delta_x) > 1e-6 ? delta_x : 1e-6);
          //take derivative of polynomial (which is coeffs[1]) to find epsi, rather than difference between first two ref points
          double epsi = -1 * atan(coeffs[1]);

          //cout <<"computed cte,epsi "<<cte<<","<<epsi<<"\n";
          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value;
          double throttle_value;
          Eigen::VectorXd state(6);
          //state << px, py, psi, v, cte, epsi;
          state << 0, 0, 0, v, cte, epsi;//in car's frame of reference, px, py, and psi are always zero
          //cout <<"calling Solve\n";
          auto vars = mpc.Solve(state, coeffs);
          steer_value = vars[0] / deg2rad(25);
          throttle_value = vars[1];
          //cur_time = std::chrono::system_clock::now();
          //elapsed_seconds = cur_time - start_time;
          //cout <<"time "<<elapsed_seconds.count()<<" : applying steer,throttle "<<steer_value<<","<<throttle_value<<"\n";

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;
          //Display the MPC predicted trajectory 
          //vector<double> mpc_x_vals;
          //vector<double> mpc_y_vals;

          vector<double> mpc_x_vals = mpc.get_x_target();
          vector<double> mpc_y_vals = mpc.get_y_target();
          vector<double> poly_x_vals;
          vector<double> poly_y_vals;
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
          // cout <<"generated points, transformed from point: ("<<cur_coord[0]<<","<<cur_coord[1]<<") and angle "<<psi<<"\n";
          //for (uint i = 0; i < raw_mpc_x_vals.size(); i++){
          //   abs_coord[0] = raw_mpc_x_vals[i];
          //   abs_coord[1] = raw_mpc_y_vals[i];
          //   //abs_coord[0] = px - i*10;
          //   //abs_coord[1] = polyeval(coeffs,abs_coord[0]);
            
          //   rel_coord = rel_coord_transform(abs_coord,cur_coord,psi);
          //   mpc_x_vals.push_back(rel_coord[0]);
          //   mpc_y_vals.push_back(rel_coord[1]);
          //   cout << "\t("<<abs_coord[0]<<","<<abs_coord[1]<<")\t("<<rel_coord[0]<<","<<rel_coord[1]<<")\n";
          // }
          //int poly_step = 3;
          //int poly_i = 0;
          //for (uint i = 0; i < mpc_x_vals.size(); i++){
          //while(poly_step * poly_i < ref_x_vals[ref_x_vals.size() - 1]){
          //  poly_x_vals.push_back(poly_i * poly_step);
          //  poly_y_vals.push_back(polyeval(coeffs,poly_x_vals[poly_i]));
          //  poly_i++;
          //}

          //plot the polynomial i'm trying to match, but with the same x values as the waypoints
          //that's too low resolution to see accuracy of polynomial.  Use 2m increments over same range
          for (double d = ref_x_vals[0]; d < ref_x_vals[ref_x_vals.size()-1]; d+=2){
            poly_x_vals.push_back(d);
            poly_y_vals.push_back(polyeval(coeffs,d));
          }
          
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;
          //msgJson["mpc_x"] = poly_x_vals;
          //msgJson["mpc_y"] = poly_y_vals;


          //msgJson["next_x"] = ref_x_vals;
          //msgJson["next_y"] = ref_y_vals;
          msgJson["next_x"] = poly_x_vals;
          msgJson["next_y"] = poly_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          
          //DONE: applied_latency set to 100 on line 15 (parameterized because I use it to figure out proper compensation)
          this_thread::sleep_for(chrono::milliseconds(applied_latency));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
