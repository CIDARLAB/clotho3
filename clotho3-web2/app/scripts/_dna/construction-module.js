"use strict";

angular.module('clotho.dna').service('Construction', ['Clotho', 'DNA', 'Digest', 'PCR', '$parse', '$q', function(Clotho, DNA, Digest, PCR, $parse, $q) {

    //testing
    var dnaModules = {};
    dnaModules.DNA = DNA;
    dnaModules.Digest = Digest;
    dnaModules.PCR = PCR;

    //reactions with templates, useable for generating new construction files
    var reactions = {
        "PCR.predict" : {
            reaction : "PCR.predict",
            input : ["", []],
            readable : "PCR"
        },
        "PCR.ligate" : {
            reaction : "PCR.ligate",
            input : [[]],
            readable : "Ligation"
        },
        "Digest.digest" : {
            reaction : "Digest.digest",
            input : ["", []],
            readable : "Restriction Digest"
        },
        "Digest.gelPurify" : {
            reaction : "Digest.gelPurify",
            input : [],
            readable : "Gel Purify"
        }
    };


    var process = function(file, returnDictionary) {
        //don't alter the construction file's dictionary
        var dict = file.dictionary;

        //wait for it to be populated in case didn't check in watch to trigger process
        if (angular.isEmpty(dict)) return returnDictionary ? [] : '';

        console.log(dict);

        //maintain references, index by key
        file.dictionaryObject = {};
        var dictObject = {};

        console.log('PROCESSING CONSTRUCTION FILE', file);

        //note - returns removed values, don't assign
        _.remove(dict, function(value) {
            return !!value.computed
        });

        //console.log(file);
        //console.log(dict);
        //console.log(_.keys(dict));


        // fixme - need to update keys, since dictionary ng-repeat requires conversion to array
        /*
        _.each(dict, function(item, oldKey) {
            if (item.key && item.key != oldKey) {
                dict[item.key] = angular.copy(item);
                delete dict[oldKey]
            }
        });
        */


        //process the dictionary, make all objects
        var dictPromises = [];
        angular.forEach(dict, function(item, index) {
            if (angular.isString(item.value)) {
                //dict[key] = {key : key, computed : false, value : value};
                angular.extend(item, {computed : false})
            }


            // if requires clotho processing, process if:
            // (1) not retrieved, (2) processed value does not equal current value
            if (!!item.Clotho) {
                if (!!item.preprocess && (!item.retrieved || item.process != item.preprocess))
                {
                    //note clientSide for now
                    var parsed = $parse(item.preprocess)(dnaModules);

                    //todo check parsed, go to clotho if undefined
                    // *** add to promises properly
                    if (angular.isUndefined(parsed)) {

                    }

                    //note - extend object so don't lose focus on proprocess changing
                    dictPromises.push = angular.extend(item, {processed : item.preprocess, retrieved : true, computed : false, value : parsed});


                    //future - dictPromises[key] = Clotho.submit(value.value);
                }
            } else {
                //todo - this is a hack to handle editing Clotho-parsed objects -- probably won't handle strings
                if (!!item.preprocess) {
                    console.log(item.value);
                    item.value = angular.isObject(item.value) ?
                        item.value :
                        JSON.parse(JSON.stringify(item.value));
                }
            }
        });


        //kickoff
        var fullDeferred = $q.defer();

        $q.all(dictPromises)
        .then(function(processedDict) {
            //angular.extend(dict, processedDict);
            //console.log(dict);
            console.log(processedDict);

            console.log(dict);

            _.each(dict, function(item) {
                dictObject[item.key] = item;
            });

            console.log(dictObject);

            var stepsChain = $q.when();

            angular.forEach(file.steps, function(step, stepIndex) {

                //run function next in chain (so dict populated)
                stepsChain.then(function() {

                    //new deferred for each step -- or return something in Clotho.run block
                    var stepDeferred = $q.defer();

                    //parse inputs
                    //future - handle promises in inputs
                    //todo - handle Clotho specifics after parse (e.g. enzyme)
                    var inputs = [];
                    angular.forEach(step.input, function(input){

                        //console.log('\n\n\n\ninput:', input);

                        //todo - better handling of empty $parse -- should still show up in computed keys

                        //future - handle 2+ layers deep
                        //for now, go one layer deep
                        if (!angular.isString(input)) {
                            //assume array or object
                            var parsedInput = (angular.isArray(input)) ? [] : {};
                            angular.forEach(input, function(child) {

                                var parsed = $parse(child)(dictObject) || {};
                                //console.log(child, parsed, dict[child]);

                                //todo - handle object scenario (don't push)
                                //todo - (later) handle defining object keys, e.g. blah[0] or blah.key

                                //for non-dictionary values, or object defined
                                //e.g. blah.value[0]
                                parsed = parsed.value || parsed;

                                parsedInput.push(parsed);
                            });
                            //console.log(parsedInput);
                            inputs.push(parsedInput);
                        } else {
                            parsedInput = $parse(input)(dictObject) || '';
                            //for non-dictionary values e.g. a boolean or number
                            parsedInput = parsedInput.value || parsedInput;
                            //console.log(parsedInput);
                            inputs.push(parsedInput);
                        }

                    });

                    //console.log(dictObject);
                    //console.log(inputs);

                    //console.log('running ' + step.reaction + ' with inputs:', inputs);

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
                    //console.log('result of '+step.reaction+' (#'+stepIndex+') with inputs:', inputs, result);
                    var newItem = {key : step.output, computed : true, stepNum: stepIndex, value : result};
                    dict.push(newItem);
                    dictObject[newItem.key] = newItem;
                    stepDeferred.resolve();

                    return stepDeferred.promise;
                });
            });

            return stepsChain;
        })
        .then(function() {

            file.dictionaryObject = dictObject;

            //console.log('final result: ', dict.final.value, dict);
            fullDeferred.resolve((!!returnDictionary) ? dict : dict.final.value);
        });

        return fullDeferred.promise;
    };

    return {
        reactions : reactions,

        process : process
    }


}]);

//note - returns array, need to maintain object references
angular.module('clotho.dna').filter('orderByProp', function() {
    return function(obj, prop) {
        return _.sortBy(obj, function (inner) {
            return inner.prop
        });
    }
});

angular.module('clotho.dna').filter('objAddKeys', function() {
    return function(obj) {
        return _.each(obj, function (value, key) {
            _.extend(value, {key : key})
        });
    }
});

angular.module('clotho.dna').filter('stepCheck', function() {
    return function(input, index) {
        return _.pick(input, function(value, key) {
            return value.computed == false || value.stepNum < index;
        });
    }
});

angular.module('clotho.dna').directive('constructionDictionaryView', [function() {
    return {
        restrict : 'EA',
        require : 'ngModel',
        scope : {
            dictionary : '=ngModel',
            editable : '=constructionEditable',
            uptostep : '=constructionUptostep',
            processto : '=constructionProcessto'
        },
        templateUrl : 'views/_dna/construction/dictionary.html',
        link : function(scope, element, attrs, ngModel) {

            scope.uptofilter = function(item) {
                return !item.computed || (item.stepNum < scope.uptostep && item.stepNum < scope.processto)
            };

            scope.addTerm = function() {
                //custom fields will be added in on process
                //todo - ensure proper form
                console.log(scope.dictionary);
                scope.dictionary.push({key : "", value :""});
            };

            /*
            //deep watch
            scope.$watch('dictionary', function(newval, oldval) {
                if (!newval || !oldval) return;
            }, true)
            */
        }
    }
}]);

angular.module('clotho.dna').directive('constructionField', ['$compile', '$filter', function($compile, $filter) {
    return {
        restrict : "EA",
        require: "ngModel",
        scope: {
            fieldName : '@constructionField',
            model : '=ngModel',
            dictionary : '=?constructionDictionary',
            options : '=constructionOptions',
            editable : '=constructionEditable',
            stepIndex : '=constructionStepindex',
            required : '@constructionRequired',
            type : '@constructionType'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {

                    var template = angular.element('<div class="control-group" ng-class="{\'error\' : hasInvalidItem || hasFutureItem}">' +
                        '<label class="control-label">{{ fieldName }}</label>' +
                        '<div class="controls">' +
                        '<constructionFieldInput></constructionFieldInput>' +
                        '</div>' +
                        '</div>');

                    template.find('constructionFieldInput').html(selectTemplate(scope.type));

                    function selectTemplate(type) {
                        switch (angular.lowercase(type)) {
                            case 'value' : {
                                return '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-disabled="!editable" ng-required="required"/>';
                            }
                            case 'boolean' : {
                                return '<button btn-checkbox class="btn" ng-class="{\'btn-success\' : !!model}"  ng-model="model" ng-disabled="!editable"><i ng-class="{\'icon-white icon-ok-circle\' : !!model, \'icon-ban-circle\' : !model}"></i></button>';
                            }
                            case 'array' : {
                                if (angular.isString(scope.model))
                                    scope.model = [scope.model];
                                return '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-list ng-change="checkInput()" ng-disabled="!editable" ng-required="required"/>';
                            }
                            default : {
                                return '<select class="span12" placeholder="{{fieldName}}" ng-model="model" ng-options="key as key for (key,val) in dict" ng-disabled="!editable" ng-change="checkInput()" ng-required="required"></select>'
                            }
                        }
                    }

                    element.html($compile(template)(scope));

                },
                post: function postLink(scope, element, attrs, ngModel) {

                    //console.log(scope.required);
                    //console.log(scope.fieldName, 'is target lenght', (scope.fieldName == "Target Length"));
                    scope.required = (scope.required != 'false') ? "true" : "false";
                    //console.log(scope.required, !!scope.required);

                    scope.$watch('dictionary', function(newval, oldval) {
                        //console.log(newval);
                        scope.dict = $filter('stepCheck')(scope.dictionary, scope.stepIndex);
                        //console.log(scope.dict);
                    });

                    scope.$watch('dict', function() {scope.checkInput()}, true);

                    scope.checkInput = function() {

                        //too early
                        if (!scope.model) return;

                        if (scope.type == 'boolean' || scope.type == 'value') return;

                        scope.hasInvalidItem = false;
                        scope.hasFutureItem = false;

                        function checkItem (model) {
                            if (!scope.dict[model])
                                scope.hasInvalidItem = true;
                            else if (scope.dict[model].stepNum >= scope.stepIndex)
                                scope.hasFutureItem = true;
                        }

                        if (scope.type == 'array') {
                            angular.forEach(scope.model, function(item) {
                                checkItem(item);
                            });
                        } else {
                            checkItem(scope.model)
                        }
                    };
                }
            }
        }
    }
}]);

angular.module('clotho.dna').directive('constructionStep', ['Construction', '$parse', '$compile', '$http', '$templateCache', '$filter', '$timeout', '$dialog', function(Construction, $parse, $compile, $http, $templateCache, $filter, $timeout, $dialog) {
    return {
        restrict : "EA",
        require: "ngModel",
        scope: {
            step : '=ngModel',
            dictionary : '=constructionDictionary',
            dictionaryObject : '=constructionDictionaryObject',
            index : '=constructionIndex',
            editable : '=constructionEditable',
            fields : '=constructionStep',
            removeStep : '&constructionRemove',
            processto : '=constructionProcessto'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {
                    
                    scope.reactions = Construction.reactions;
                    scope.reactionsArray = _.toArray(scope.reactions);

                    //note - can't recompile step or lose focus every keystroke dictionary model changes if force recompile in dictionary watch, so pull out into process step
                    scope.processStep = function processStep() {
                        //scope.dict = $filter('stepCheck')(scope.dictionary, scope.index);
                        scope.dict = scope.dictionaryObject;
                    };

                    scope.compileStep = function compileStep() {

                        scope.processStep();

                        var template = angular.element('<div class="constructionStep clearfix">'+
                            '<div class="arrow"></div>' +
                            '<div class="reaction" ng-class="{\'errorBackground\' : !reactions[step.reaction]}">' +
                            '<i ng-class="{\'icon-resize-full\' : !showReaction, \'icon-resize-small\' : showReaction}" style="cursor : pointer;" ng-click="toggleReaction()"></i> ' +
                            '{{ reactions[step.reaction].readable || step.reaction }}'+
                            //'<input type="text" class="reactionName" typeahead="r.reaction as r.readable for r in reactionsArray" ng-model="step.reaction" style="margin-bottom: 0px; padding: 3px 5px;" />' +
                            '<i ng-if="editable" class="pull-right icon-remove" ng-click="confirmRemoveStep()" style="cursor : pointer; margin-top: 3px"></i> ' +
                            '</div>'+
                            '<div ng-show="showReaction">' +
                            //'{{fields}}' +
                            '<div class="fields" ng-style="{cursor : (editable ? \'move\' : \'inherit\')}">' +
                            '<stepFields class="row-fluid clearfix"></stepFields>' +
                            //'{{ step }}'+
                            '</div>' +
                            '<div class="output" ng-class="{\'errorBackground\' : !dictionaryObject[step.output], \'warningBackground\' : !dictionaryObject[step.output].value}" tooltip="{{ index < processto ? (dictionaryObject[step.output].value | stringEnds) : \'<unprocessed>\' }}">Output: <code contenteditable="{{ editable }}" ng-model="step.output" style="display: inline-block; padding: 0 4px;"></code></div>' +
                            '</div>');


                        function generateConstructionField (type, name, model) {

                            return '<div class="span4" construction-field="'+name+'" construction-type="'+type+'" construction-editable="editable" construction-stepindex="index" construction-dictionary="dict" ng-model="'+model+'" ></div>';
                        }

                        //if pass in scope.fields, assume want to overrun templates
                        if (angular.isDefined(scope.fields)) {
                            angular.forEach(scope.fields, function(field) {
                                field.model = 'step.' + field.model;

                                var html = generateConstructionField(field.type, field.name, field.model);
                                template.find('stepFields').append(html)
                            });
                        }
                        else {
                            //try to get the template for the reaction
                            $http.get('views/_dna/construction/steps/' + scope.step.reaction + '.html', {cache : $templateCache})
                                .then(function(data) {
                                    //works - great
                                    //console.log('got template for ' + scope.step.reaction + '\n', data.data)
                                    //console.log(data.data);
                                    //console.log(scope);
                                    var fieldsTemplate = $compile(data.data)(scope);
                                    template.find('stepFields').append(fieldsTemplate);

                                }, function(data) {
                                    //didn't  work - try to parse input fields
                                    console.log('template not found');

                                    angular.forEach(scope.step.input, function(input, index) {
                                        var html;

                                        if (angular.isArray(input)) {
                                            html = generateConstructionField('array', '', 'step.input['+index+']');
                                        }
                                        else if (input == "true" || input == "false") {
                                            html = generateConstructionField('boolean', '', 'step.input['+index+']');
                                        }
                                        else {
                                            html = generateConstructionField('value', '', 'step.input['+index+']');
                                        }

                                        template.find('stepFields').append(html)
                                    });
                                });
                        }

                        element.html($compile(template)(scope))
                    };

                    scope.compileStep();

                },
                post: function postLink(scope, element, attrs, ngModel) {

                    //future - alternatively, filter to disable select options
                        // http://jsfiddle.net/alalonde/dZDLg/9/

                    //watch dictionary so runs after process
                    /*scope.$watch('dictionary', function (newval, oldval) {
                        scope.processStep();
                    });*/
                    scope.$watch('dictionaryObject', function (newval, oldval) {
                        scope.processStep();
                    });

                    scope.showReaction = true;
                    scope.toggleReaction = function() {
                        scope.showReaction = !scope.showReaction;
                    };

                    scope.confirmRemoveStep = function() {
                        var d = $dialog.messageBox('Confirm Remove', 'Are you sure you want to remove this reaction?', [{label: "Cancel", cssClass: "", result: false}, {label: "Remove", cssClass: "btn-danger", result: true}]);
                        d.open().then(function (result) {
                            if (result) {
                                scope.removeStep();
                            }
                        });
                    }
                }
            }
        }
    }
}]);

angular.module('clotho.dna').directive('constructionFull', ['Construction', function(Construction) {
    return {
        restrict : 'A',
        require : 'ngModel',
        replace : true,
        scope : {
            file : '=ngModel',
            editable : '=constructionEditable',
            onChange : '&?constructionOnchange',
            hideopts : '=?constructionHideopts'
        },
        templateUrl : "views/_dna/construction/full.html",
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, controllers) {

                },
                post: function postLink(scope, element, attrs, controllers) {

                    var hideOptsDef = {
                        dictionary : false,
                        steps : false,
                        addStep : false,
                        uptostep : 99,
                        processto : 99
                    };
                    scope.hideopts = angular.extend(hideOptsDef, scope.hideopts);

                    scope.$watch('file', function(newval, oldval) {
                        if (!newval) return;
                        Construction.process(newval, true).then(function(result) {
                            console.log('finalResult', result);

                            //need to update dictionary since refernece not maintained in processing, forces another digest
                            scope.file.dictionary = result;
                            //scope.product = scope.file.dictionary.final.value;

                            scope.onChange({file : scope.file});
                        });
                    }, true);

                    scope.reactions = Construction.reactions;
                    scope.reactionsArray = _.toArray(scope.reactions);

                    scope.newStep = {
                        reaction : "",
                        input : [],
                        output : ""
                    };
                    scope.stepValid = function() {
                        return !!Construction.reactions[scope.newStep.reaction]
                    };
                    scope.addStep = function() {
                        if (scope.stepValid(scope.newStep)) {
                            scope.newStep.input = Construction.reactions[scope.newStep.reaction].input;
                            scope.file.steps.push(scope.newStep);
                            scope.newStep = {};
                            scope.hideopts.uptostep = scope.file.steps.length;
                        }
                    }

                }
            }
        }
    }
}]);