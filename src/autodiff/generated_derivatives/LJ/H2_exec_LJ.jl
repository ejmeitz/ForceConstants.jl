function H2_exec_LJ(input_variables)
# input_variables->begin
        #= C:\Users\ejmei\.julia\packages\FastDifferentiation\mftkg\src\CodeGeneration.jl:192 =#
        #= C:\Users\ejmei\.julia\packages\FastDifferentiation\mftkg\src\CodeGeneration.jl:192 =# @inbounds begin
                #= C:\Users\ejmei\.julia\packages\FastDifferentiation\mftkg\src\CodeGeneration.jl:193 =#
                begin
                    result_element_type = promote_type(Float64, eltype(input_variables))
                    result = Array{result_element_type}(undef, (3, 3))
                    var"##965" = input_variables[1] ^ 2
                    var"##966" = input_variables[2] ^ 2
                    var"##964" = var"##965" + var"##966"
                    var"##967" = input_variables[3] ^ 2
                    var"##963" = var"##964" + var"##967"
                    var"##962" = sqrt(var"##963")
                    var"##961" = 2var"##962"
                    var"##960" = 1 / var"##961"
                    var"##973" = 3.4 / var"##962"
                    var"##972" = var"##973" / var"##962"
                    var"##971" = -var"##972"
                    var"##977" = var"##973" ^ 5
                    var"##976" = var"##977" * var"##977"
                    var"##975" = 69.22656var"##976"
                    var"##980" = var"##973" ^ 4
                    var"##984" = var"##973" ^ 6
                    var"##983" = var"##984" - 1
                    var"##982" = 0.96148var"##983"
                    var"##985" = 0.96148var"##984"
                    var"##981" = var"##982" + var"##985"
                    var"##979" = var"##980" * var"##981"
                    var"##978" = 30var"##979"
                    var"##974" = var"##975" + var"##978"
                    var"##970" = var"##971" * var"##974"
                    var"##969" = var"##970" * var"##971"
                    var"##989" = 6var"##977"
                    var"##988" = var"##981" * var"##989"
                    var"##987" = -var"##988"
                    var"##992" = 1 / var"##962"
                    var"##991" = var"##992" * var"##971"
                    var"##994" = var"##972" / var"##962"
                    var"##993" = -var"##994"
                    var"##990" = var"##991" + var"##993"
                    var"##986" = var"##987" * var"##990"
                    var"##968" = var"##969" + var"##986"
                    var"##959" = var"##960" * var"##968"
                    var"##997" = var"##988" * var"##971"
                    var"##999" = var"##960" / var"##961"
                    var"##998" = -var"##999"
                    var"##996" = var"##997" * var"##998"
                    var"##995" = 2var"##996"
                    var"##958" = var"##959" + var"##995"
                    var"##957" = var"##958" * var"##960"
                    var"##1001" = input_variables[1] * input_variables[1]
                    var"##1000" = 4var"##1001"
                    var"##956" = var"##957" * var"##1000"
                    var"##1003" = var"##997" * var"##960"
                    var"##1002" = 2var"##1003"
                    var"##955" = var"##956" + var"##1002"
                    result[CartesianIndex(1, 1)] = var"##955"
                    var"##1007" = 2 * input_variables[2]
                    var"##1006" = var"##1007" * var"##958"
                    var"##1005" = var"##1006" * var"##960"
                    var"##1008" = 2 * input_variables[1]
                    var"##1004" = var"##1005" * var"##1008"
                    result[CartesianIndex(2, 1)] = var"##1004"
                    var"##1012" = 2 * input_variables[3]
                    var"##1011" = var"##1012" * var"##958"
                    var"##1010" = var"##1011" * var"##960"
                    var"##1009" = var"##1010" * var"##1008"
                    result[CartesianIndex(3, 1)] = var"##1009"
                    var"##1016" = input_variables[2] * input_variables[1]
                    var"##1015" = 4var"##1016"
                    var"##1014" = var"##1015" * var"##958"
                    var"##1013" = var"##1014" * var"##960"
                    result[CartesianIndex(1, 2)] = var"##1013"
                    var"##1020" = input_variables[2] * input_variables[2]
                    var"##1019" = 4var"##1020"
                    var"##1018" = var"##957" * var"##1019"
                    var"##1017" = var"##1018" + var"##1002"
                    result[CartesianIndex(2, 2)] = var"##1017"
                    var"##1021" = var"##1010" * var"##1007"
                    result[CartesianIndex(3, 2)] = var"##1021"
                    var"##1025" = input_variables[3] * input_variables[1]
                    var"##1024" = 4var"##1025"
                    var"##1023" = var"##1024" * var"##958"
                    var"##1022" = var"##1023" * var"##960"
                    result[CartesianIndex(1, 3)] = var"##1022"
                    var"##1026" = var"##1005" * var"##1012"
                    result[CartesianIndex(2, 3)] = var"##1026"
                    var"##1030" = input_variables[3] * input_variables[3]
                    var"##1029" = 4var"##1030"
                    var"##1028" = var"##957" * var"##1029"
                    var"##1027" = var"##1028" + var"##1002"
                    result[CartesianIndex(3, 3)] = var"##1027"
                    return result
                end
            end
    end