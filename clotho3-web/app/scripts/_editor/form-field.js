angular.module('clotho.editor')
/**
 * @name formField
 *
 * @description Wrapper for form elements, adding bootstrap classes automatically. Also adds Labels and Help blocks. Must be recompiled so don't change internal contents and expect re-render.
 *
 * @note Known problem: Updates to form controls within this directive will not an undefined model. Model should at least be declared as empty object.
 *
 * @usage Use in form on elements: input, textarea, select, etc. Not button. Should have one and only one direct child (the field element). That child may have it's own children (e.g. options in a select).
 *
 * @example
 <form-field name="Long Input" help="Enter a short biography about yourself">
 <textarea rows="3" id="exampleTextarea" ng-model="model.bio" placeholder="Write a short biography"></textarea>
 </form-field>

 Produces:

 <div ng-form="" class="form-group" name="Long Input" help="Enter a short biography about yourself">
 <label for="exampleTextarea">Long Input</label>
 <textarea rows="3" id="exampleTextarea" ng-model="model.bio" placeholder="Write a short biography" class="form-control"></textarea>
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
			transclude: true,
			require: '^form',
			controller: function ($scope, $element, $attrs) {

			},
			link : function linkFunction(scope, element, attrs, formCtrl) {

				var passedName = attrs.name,
					passedHelp = attrs.help,
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
					wrapper.append(label);
					element.append(wrapper);
				}

				//add label if pass name
				else if (passedName) {
					element.prepend(label);
				}



				//append help
				if (passedHelp) {
					element.append('<p class="help-block">'+passedHelp+'</p>');
				}

				//watch for error
				scope.$watch(function() {
					return formCtrl.$invalid && formCtrl.$error;
				}, function (isValid, lastVal) {
					element.toggleClass('has-error', isValid);
				});
			}
		};
	});