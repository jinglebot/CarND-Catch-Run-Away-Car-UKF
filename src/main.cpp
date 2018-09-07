#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

int main()
{
  uWS::Hub h;

  // Create a UKF instance
  UKF ukf;

  double target_x = 0.0;
  double target_y = 0.0;
  double target_v = 0.0;

  h.onMessage([&ukf,&target_x,&target_y,&target_v](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event

    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {

      auto s = hasData(std::string(data));
      if (s != "") {
      	      	
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object

          double hunter_x = std::stod(j[1]["hunter_x"].get<std::string>());
          double hunter_y = std::stod(j[1]["hunter_y"].get<std::string>());
          double hunter_heading = std::stod(j[1]["hunter_heading"].get<std::string>());
          
          string lidar_measurment = j[1]["lidar_measurement"];
          
          MeasurementPackage meas_package_L;
          istringstream iss_L(lidar_measurment);
      	  long long timestamp_L;

      	  // reads first element from the current line
      	  string sensor_type_L;
      	  iss_L >> sensor_type_L;

      	  // read measurements at this timestamp
      	  meas_package_L.sensor_type_ = MeasurementPackage::LASER;
          meas_package_L.raw_measurements_ = VectorXd(2);
          float px;
      	  float py;
          iss_L >> px;
          iss_L >> py;
          meas_package_L.raw_measurements_ << px, py;
          iss_L >> timestamp_L;
          meas_package_L.timestamp_ = timestamp_L;
          
      	  ukf.ProcessMeasurement(meas_package_L);
         
      	  string radar_measurment = j[1]["radar_measurement"];
          
          MeasurementPackage meas_package_R;
          istringstream iss_R(radar_measurment);
      	  long long timestamp_R;

      	  // reads first element from the current line
      	  string sensor_type_R;
      	  iss_R >> sensor_type_R;

      	  // read measurements at this timestamp
      	  meas_package_R.sensor_type_ = MeasurementPackage::RADAR;
          meas_package_R.raw_measurements_ = VectorXd(3);
          float ro;
      	  float theta;
      	  float ro_dot;
          iss_R >> ro;
          iss_R >> theta;
          iss_R >> ro_dot;
          meas_package_R.raw_measurements_ << ro,theta, ro_dot;
          iss_R >> timestamp_R;
          meas_package_R.timestamp_ = timestamp_R;
          
      	  ukf.ProcessMeasurement(meas_package_R);
        
          target_x = ukf.x_[0];
          target_y = ukf.x_[1];
          target_v = ukf.x_[2];
          

          double future_x;
          double future_y;   
          double distance_difference;
          ukf.iteration++;

          // after 250 iterations it is assumed the error covariance is fully minimized
          if (ukf.iteration > 250) {

            double delta_t = 0;
            double current_difference = 0;        
            
            double previous_x = target_x;
            double previous_y = target_y;

            double previous_distance_difference = sqrt((target_y - hunter_y)*(target_y - hunter_y) + (target_x - hunter_x)*(target_x - hunter_x));
            double previous_difference = fabs(previous_distance_difference - fabs(target_v * delta_t));

            do {
              delta_t += 0.02;
              // distance from target current state to the future state 
              if (fabs(ukf.x_[4]) > 0.0001) {
                future_x = target_x + (target_v / ukf.x_[4]) * (sin(ukf.x_[3] + ukf.x_[4] * delta_t) - sin(ukf.x_[3]));
                future_y = target_y + (target_v / ukf.x_[4]) * (cos(ukf.x_[3]) - cos(ukf.x_[3] + ukf.x_[4] * delta_t));
              } else {
                future_x = target_x + target_v * delta_t * cos(ukf.x_[3]);
                future_y = target_y + target_v * delta_t * sin(ukf.x_[3]);
              }
              // ukf_future.Prediction(delta_t); tried it, didn't work
              // distance from intercept position to hunter position
              distance_difference = sqrt((future_y - hunter_y)*(future_y - hunter_y) + (future_x - hunter_x)*(future_x - hunter_x));
              // double target_difference = sqrt((future_y - target_y)*(future_y - target_y) + (future_x - target_x)*(future_x - target_x)); tried it, didn't look right
              // you want the distance the target travelled to intercept == that of hunter or their difference at a minimum
              current_difference = fabs(distance_difference - fabs(target_v * delta_t));
              // if the distance to the intercept point from the hunter == that from the target
              if (current_difference < 0.09) {
                break;
                std::cout << "current_difference hunter to future == target to future" << std::endl;    
              }
              // if the current difference error is bigger, don't change the values, use the previous
              if (current_difference > previous_difference) {
                // hunter stays on the same path using previous values
                future_x = previous_x;
                future_y = previous_y;
                std::cout << "previous values used" << std::endl;            
              }
              // update previous values
              previous_x = future_x;
              previous_y = future_y;
              // maximum time duration to calculate values == 5 secs 
            } while (delta_t < 3);

            std::cout << "delta_t: " << delta_t << std::endl;            
            std::cout << "difference: " << current_difference << std::endl;
            std::cout << "**************************" << std::endl;            
          }

          double heading_to_target = atan2(future_y - hunter_y, future_x - hunter_x);
          ukf.Normalize_angle(heading_to_target);
          //turn towards the target
          double heading_difference = heading_to_target - hunter_heading;
          ukf.Normalize_angle(heading_difference);

          json msgJson;
          // use multipliers to change a little of the heading and distance per iteration
          msgJson["turn"] = heading_difference * 0.05;
          msgJson["dist"] = distance_difference * 0.4; 
          auto msg = "42[\"move_hunter\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
	  
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }

  });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}









































































