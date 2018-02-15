#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

class MeasurementPackage {
public:
  long timestamp_;

  enum SensorType{
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;

  // Keep a record of last measured position for the logging purpose
  double measured_px_;
  double measured_py_;

  void UpdatePosition() {
    if (sensor_type_ == MeasurementPackage::RADAR) {
      double ro = raw_measurements_(0);
      double theta = raw_measurements_(1);
      measured_px_ = ro * cos(theta);
      measured_py_ = ro * sin(theta);
    }
    else if (sensor_type_ == MeasurementPackage::LASER) {
      measured_px_ = raw_measurements_(0);
      measured_py_ = raw_measurements_(1);
    }
  }
};

#endif /* MEASUREMENT_PACKAGE_H_ */
