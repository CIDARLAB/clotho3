"use strict";

Application.Dna.service('Construction', ['Clotho', 'DNA', 'Digest', 'PCR', '$parse', '$q', function(Clotho, DNA, Digest, PCR, $parse, $q) {

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

        console.log('PROCESSING CONSTRUCTION FILE', file);

        dict = _.pick(dict, function(value, key) {
            return !value.computed
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
        var dictPromises = {};
        angular.forEach(dict, function(value, key) {
            if (angular.isString(value)) {
                dict[key] = {key : key, computed : false, value : value};
            }

            if (!!value.Clotho && !(value.retrieved === true)) {
                //testing clientSide for now
                var parsed = $parse(value.value)(dnaModules);
                dictPromises[key] = {key: key, Clotho : true, retrieved : true, computed : false, value : parsed};

                //future - dictPromises[key] = Clotho.submit(value.value);
            }
        });


        //kickoff
        var fullDeferred = $q.defer();

        $q.all(dictPromises)
        .then(function(processedDict) {
            angular.extend(dict, processedDict);

            //console.log(dict);

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

                                var parsed = $parse(child)(dict) || {};
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
                            parsedInput = $parse(input)(dict) || '';
                            //for non-dictionary values e.g. a boolean or number
                            parsedInput = parsedInput.value || parsedInput;
                            //console.log(parsedInput);
                            inputs.push(parsedInput);
                        }
                    });

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
                    dict[step.output] = {key : step.output, computed : true, stepNum: stepIndex, value : result};
                    stepDeferred.resolve();

                    return stepDeferred.promise;
                });
            });

            return stepsChain;
        })
        .then(function() {

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

//todo - returns array, need to maintain object references
Application.Dna.filter('orderByProp', function() {
    return function(obj, prop) {

        return _.sortBy(obj, function (inner) {
            return inner.prop
        });
    }
});

Application.Dna.filter('objAddKeys', function() {
    return function(obj) {
        return _.each(obj, function (value, key) {
            _.extend(value, {key : key})
        });
    }
});

Application.Dna.filter('stepCheck', function() {
    return function(input, index) {
        return _.pick(input, function(value, key) {
            return value.computed == false || value.stepNum < index;
        });
    }
});

Application.Dna.directive('constructionDictionaryView', [function() {
    return {
        restrict : 'EA',
        require : 'ngModel',
        scope : {
            dictionary : '=ngModel',
            editable : '=constructionEditable',
            uptostep : '=constructionUptostep'
        },
        template : '<div style="overflow: scroll">' +
            //'<pre class="pre-scrollable" ng-bind="dictionary | json"></pre>' +
            '<table class="table table-bordered table-condensed constructionDictionary">' +
            '<thead><tr><th>Key</th> <th>Value</th></tr></thead>' +
            '<tbody>' +
                '<tr ng-repeat="item in dictionary | orderByProp:\'stepNum\' | filter:uptofilter" ng-class="{\'computed\' : item.computed}">' +
                '<td><code contenteditable="{{ !!editable && !item.computed}}" ng-model="item.key" tooltip="{{ !!editable && item.computed ? \'Edit computed keys in the workflow\' : \'\' }}"></code></td>' +
                //'<td>{{item}}</td>' +
                '<td contenteditable="{{ !!editable && !item.computed && !item.retrieved}}" tooltip="{{ !!editable && item.computed ? \'You cannot edit a dynamically generated value\' : \'\' }}" tooltip-placement="mouse" ng-model="item.value" style="overflow: hidden; text-overflow: ellipsis;"></td>' +
                '</tr>' +
                '<tr ng-if="!!editable"><td colspan="2"><button class="btn" ng-click="addTerm()">Add new Key</button></td></tr>' +
            '</tbody>' +
            '</table>' +
            '</div>',
        link : function(scope, element, attrs, ngModel) {

            scope.uptofilter = function(item) {
                return !item.computed || item.stepNum < scope.uptostep
            };

            scope.addTerm = function() {
                //custom fields will be added in on process
                //todo -- add ability to clotho parse
                angular.extend(scope.dictionary, {"" : ""});
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

Application.Dna.directive('constructionField', ['$compile', '$filter', function($compile, $filter) {
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

Application.Dna.directive('constructionStep', ['Construction', '$parse', '$compile', '$http', '$templateCache', '$filter', '$timeout', '$dialog', function(Construction, $parse, $compile, $http, $templateCache, $filter, $timeout, $dialog) {
    return {
        restrict : "EA",
        require: "ngModel",
        scope: {
            step : '=ngModel',
            dictionary : '=constructionDictionary',
            index : '=constructionIndex',
            editable : '=constructionEditable',
            fields : '=constructionStep',
            removeStep : '&constructionRemove'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {

                    scope.reactions = Construction.reactions;
                    scope.reactionsArray = _.toArray(scope.reactions);

                    //note - can't recompile step or lose focus every keystroke dictionary model changes if force recompile in dictionary watch, so pull out into process step
                    scope.processStep = function processStep() {
                        //scope.dict = $filter('stepCheck')(scope.dictionary, scope.index);
                        scope.dict = scope.dictionary;
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
                            '<div class="output" ng-class="{\'errorBackground\' : !dictionary[step.output], \'warningBackground\' : !dictionary[step.output].value}" tooltip="{{ dictionary[step.output].value | stringEnds }}">Output: <code contenteditable="{{ editable }}" ng-model="step.output" style="display: inline-block; padding: 0 4px;"></code></div>' +
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
                            $http.get('/dna/construction/steps/' + scope.step.reaction + '.html', {cache : $templateCache})
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
                    scope.$watch('dictionary', function (newval, oldval) {
                        scope.processStep();

                        /*testing
                        _.each(newval, function(item){
                           console.log(scope.index, item.key, item.computed, item.stepNum);
                        });
                        */
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

Application.Dna.directive('constructionFull', ['Construction', function(Construction) {
    return {
        restrict : 'A',
        require : 'ngModel',
        replace : true,
        scope : {
            file : '=ngModel',
            editable : '=constructionEditable',
            onChange : '&?constructionOnchange',
            hideopts : '=constructionHideopts'
        },
        templateUrl : "/dna/construction/full.html",
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, controllers) {



                },
                post: function postLink(scope, element, attrs, controllers) {

                    var hideOptsDef = {
                        dictionary : false,
                        steps : false,
                        addStep : false,
                        uptostep : 99
                    };
                    scope.hideopts = angular.extend(hideOptsDef, scope.hideopts);

                    scope.$watch('file', function(newval, oldval) {
                        if (!newval) return;
                        Construction.process(newval, true).then(function(result) {
                            console.log('finalResult', result);

                            //need to update dictionary since refernece not maintained in processing, forces another digest
                            scope.file.dictionary = result;
                            scope.product = scope.file.dictionary.final.value;

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
}])