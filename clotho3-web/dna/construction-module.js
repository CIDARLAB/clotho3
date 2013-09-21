"use strict";

Application.Dna.service('Construction', ['Clotho', 'DNA', 'Digest', 'PCR', '$parse', '$q', function(Clotho, DNA, Digest, PCR, $parse, $q) {

    //testing
    var dnaModules = {};
    dnaModules.DNA = DNA;
    dnaModules.Digest = Digest;
    dnaModules.PCR = PCR;

    var process = function(file) {
        //don't alter the construction file's dictionary
        var dict = angular.copy(file.dictionary);

        //process the dictionary
        var dictPromises = {};
        angular.forEach(dict, function(value, key) {
            if (!!value.Clotho) {
                //dictPromises[key] = Clotho.submit(value.value);

                //testing clientSide for now
                dictPromises[key] = $parse(value.value)(dnaModules);
            }
        });


        var stepsChain = $q.defer();
        angular.forEach(file.steps, function(step) {

            //run function next in chain (so dict populated)
            stepsChain.promise.then(function() {

                //new deferred for each step -- or return something in Clotho.run block
                var stepDeferred = $q.defer();

                //parse inputs
                //future - handle promises in inputs
                //todo - handle Clotho specifics after parse (e.g. enzyme)
                var inputs = [];
                angular.forEach(step.inputs, function(input){
                    //future - handle 2+ layers deep
                    //for now, go one layer deep
                    if (!angular.isString(input)) {
                        //assume array or object
                        var parsedInput = (angular.isArray(input)) ? [] : {};
                        angular.forEach(input, function(child) {
                            //todo - handle object scenario
                            parsedInput.push($parse(child)(dict));
                        });
                        inputs.push(parsedInput);
                    } else {
                        inputs.push($parse(input)(dict));
                    }
                });

                console.log('running ' + step.reaction + ' with inputs:', inputs);

                //future -  run in clotho
                /*
                Clotho.run(step.reaction, inputs).then(function(result) {
                    dict[step.output] = result;
                    console.log(dict);
                    stepDeferred.resolve()
                });
                */

                //note - for now running client side
                //use apply to match format on server
                var result = $parse(step.reaction)(dnaModules).apply(null, inputs);
                console.log('result of '+step.reaction+':', result);
                dict[step.output] = result;

                return stepDeferred.promise;
            });
        });

        //kickoff
        $q.all(dictPromises).then(function(processedDict) {
            angular.extend(dict, processedDict);
            stepsChain.resolve();
        });

        //at the end
        return stepsChain.promise.then(function() {
            console.log('final: ' + dict.final, dict);
            return dict.final
        });


    };

    return {
        process : process
    }


}]);