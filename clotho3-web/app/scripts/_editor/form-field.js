angular.module('clotho.editor')
/**
 * @name formField
 *
 * @description Wrapper for form elements, adding bootstrap classes automatically. Specifically handles according to element type, wrapping appropriately. Also adds Labels and Help blocks. Does not compile internal contents.
 *
 * @note The element declared as child of form-field is transcluded, ng-model delcarations will run into prototypical inheritance weirdness (i.e. if you don't declare in the form of "object.field", define as "$parent.object"). Updates to form controls within this directive will not an undefined model. Model should at least be declared as empty object.
 *
 * @usage Use in form on elements: input, textarea, select, etc. Not button. Should have one and only one direct child (the field element). That child may have it's own children (e.g. options in a select). If pass removable="fieldName" then will delete key sharable.fieldName. Can pass hide-label to hide label and just show the element and help, but will be overriden for input types that surround the input with the label (e.g. checkbox, radio)
 *
 * @example
 <form-field name="Long Input" help="Enter a short biography about yourself" hide-label="false" removable="true">
 <textarea rows="3" id="exampleTextarea" ng-model="model.bio" placeholder="Write a short biography"></textarea>
 </form-field>

 Produces (some removeField classes stripped):

 <div ng-form="" class="form-group" name="Long Input" help="Enter a short biography about yourself">
 <label for="exampleTextarea">Long Input</label>
 <textarea rows="3" id="exampleTextarea" ng-model="model.bio" placeholder="Write a short biography" class="form-control"></textarea>
 <button ng-click="removeField($index)"></button>
 <p class="help-block">Enter a short biography about yourself</p>
 </div>

 */
	.directive('formField', function () {

		var template = '<div class="form-group" ng-form ng-transclude>' +
			'</div>';

		return {
			restrict: 'E',
			template: template,
			replace: true,
			transclude: true, //only want element contents
			require: ['^form', '^?clothoEditor'],
			controller: function ($scope, $element, $attrs) {

			},
			link : function linkFunction(scope, element, attrs, parentCtrls) {

				var formCtrl = parentCtrls[0],
					editorCtrl = parentCtrls[1];

				var passedName = attrs.name,
					passedHelp = attrs.help,
					removable = attrs.removableField,
					hideLabel = attrs.hideLabel,
					childElement = element.children();

				if (childElement.length !== 1) {
					throw 'You must include the form input element as the single child of this directive, instead got ' + element.html() + 'generating:';
				}

				var elemId = childElement.attr('id');
				var elemTag = childElement.prop('tagName');
				var elemType = childElement.attr('type');

				//regex for element types that are special
				var regWrapElementInType = /checkbox|radio/gi;
				var regNoFormCtrlClass = /file|checkbox|radio/gi;

				//if didn't set element id, then generate and set ourselves
				if (!elemId) {
					elemId = 'input' + (Math.floor(Math.random()*10000000)).toString();
					childElement.attr('id', elemId);
				}

				//child class setup
				if (!regNoFormCtrlClass.test(elemType)){
					childElement.addClass('form-control');
				}

				//label
				var label = angular.element('<label class="control-label" for="'+elemId+'">'+
					( elemType === 'radio' ? childElement.val() : passedName) +
					'</label>');

				//wrap with bootstrap class if appropriate
				if (regWrapElementInType.test(elemType)) {
					label.prepend(childElement);
					var wrapper = angular.element('<div class="'+elemType+'"></div>');
					element.append(wrapper);
					hideLabel = false;
				}

				/*
				//note - not supporting field removal for now
				scope.removeField = angular.isDefined(editorCtrl) ? editorCtrl.removeField : angular.noop;
				if (removable) {
					var wrapper = angular.element('<div class="input-group"></div>');
					var removeButton = angular.element('<span class="input-group-btn"><button class="btn btn-danger" type="button" ng-click="removeField(\''+removable+'\')"><span class="glyphicon glyphicon-trash"></span></button></span>');
					$compile(removeButton)(scope);
					wrapper.prepend(element.contents());
					wrapper.append(removeButton);
					element.append(wrapper);
				}
			 */

				//add label if pass name
				if (passedName && !hideLabel) {
					element.prepend(label);
				}

				//append help
				if (passedHelp) {
					element.append('<p class="help-block">'+passedHelp+'</p>');
				}

				//watch for error
				scope.$watch(function() {
					return formCtrl.$invalid;
				}, function (isInvalid, lastVal) {
					element.toggleClass('has-error', isInvalid);
				});
			}
		};
	});