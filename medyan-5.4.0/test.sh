#Run unit tests
echo "unit tests"
./build/medyan test

echo "2filaments"
DIR=2filaments
mkdir -p ./examples/${DIR}/out
./build/medyan -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "50filaments_motor_linker"
DIR=50filaments_motor_linker
mkdir -p ./examples/${DIR}/out
./build/medyan -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "actin only"
DIR=actin_only
mkdir -p ./examples/${DIR}/out
./build/medyan -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "branch_actin"
DIR=branch_actin
mkdir -p ./examples/${DIR}/out
./build/medyan -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "membrane_vesicle"
DIR=membrane_vesicle
mkdir -p ./examples/${DIR}/out
./build/medyan -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "nucleation_actin"
DIR=nucleation_actin
mkdir -p ./examples/${DIR}/out
./build/medyan -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out