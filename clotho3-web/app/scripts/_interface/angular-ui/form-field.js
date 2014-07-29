angular.module('clotho.interface')
/**
 * @name formField
 *
 * @description Wrapper for form elements, adding bootstrap classes automatically. Specifically handles according to element type, wrapping appropriately. Also adds Labels and Help blocks. Does not compile internal contents.
 *
 * @note Models will not propagate if they are primitives. Use objects / array with a dot (e.g. my.model)
 *
 * The element declared as child of form-field is transcluded, ng-model delcarations may run into prototypical inheritance weirdness [transcluded elements SHUOLD NOT have isolate scope, sibling to directive, so shuold not be affected by directive in prototypical inheritance] (i.e. if you don't declare in the form of "object.field", define as "$parent.object"). Updates to form controls within this directive will not update an undefined model. Model should at least be declared as empty object / array.
 *
 * @usage Use in form on elements: input, textarea, select, etc. Not button. Should have one and only one direct child (the field element). That child may have it's own children (e.g. options in a select). If pass removable="fieldName" then will delete key sharable.fieldName. Can pass hide-label to hide label and just show the element and help, but will be overriden for input types that surround the input with the label (e.g. checkbox, radio)
 *
 * @attr name
 * @attr removable
 * @attr help
 * @attr noStyling
 *
 * @example
 <form-field name="Long Input" help="Enter a short biography about yourself"  removable="true" horizontal="true">
 <textarea rows="3" id="exampleTextarea" ng-model="model.bio" placeholder="Write a short biography"></textarea>
 </form-field>

 Produces (some removeField classes stripped):

 <ng-form class="form-group" name="Long Input" help="Enter a short biography about yourself">
	 <label for="exampleTextarea">Long Input</label>
	 <textarea rows="3" id="exampleTextarea" ng-model="model.bio" placeholder="Write a short biography" class="form-control"></textarea>
	 <button ng-click="removeField($index)"></button>
	 <p class="help-block">Enter a short biography about yourself</p>
 </ng-form>

 */
	.directive('formField', function ($parse) {

		var template = '<ng-form class="form-group" ng-transclude></ng-form>';

		return {
			restrict: 'E',
			template: template,
			replace: true,
			transclude: true, //only want element contents
			require: ['^form'],
			controller: function ($scope, $element, $attrs) {

			},
			link : function linkFunction(scope, element, attrs, parentCtrls, transclude) {

				var formCtrl = parentCtrls[0];

				var passedName = attrs.name,
					passedHelp = attrs.help,
					horizontal = $parse(attrs.horizontal)(scope),
					removable = attrs.removableField,
					labelAdded = false,
					childElement = element.children();

				if (childElement.length !== 1) {
					throw 'You must include the form input element as the single child of this directive, instead got ' + element.html() + 'generating:';
				}

				var elemId = childElement.attr('id');
				var elemTag = childElement.prop('tagName');
				var elemType = childElement.attr('type');

				//regex for element types that are special
				var regWrapElementInType = /checkbox|radio|button/gi;
				var regNoFormCtrlClass = /file|checkbox|radio/gi;

				//if didn't set element id, then generate and set ourselves
				if (!elemId) {
					elemId = 'input' + (Math.floor(Math.random()*10000000)).toString();
					childElement.attr('id', elemId);
				}

				//add form-control class to child UNLESS bootstrap says no, or noStyling attr is set
				if (! ( regNoFormCtrlClass.test(elemType) || angular.isDefined(attrs.noStyling) )) {
					console.log('adding form control class,',  angular.isDefined(attrs.noStyling), element);
					childElement.addClass('form-control');
				}

				//create label
				var label = angular.element('<label class="control-label" for="'+elemId+'">'+
					( elemType === 'radio' ? childElement.val() : passedName) +
					'</label>');

				//wrap with bootstrap class if appropriate (e.g. checkbox, radio)
				if (regWrapElementInType.test(elemType)) {
					var wrapper = angular.element('<div class="'+elemType+'"></div>');
					/*
					//for normal bootstrap way of handling, but leaves nothing in left column as label...
					label.prepend(childElement);
					childElement = wrapper.append(label);
					labelAdded = true;
					*/
					wrapper.append(childElement);
					childElement = wrapper;
				}

				//if horizontal, add label to side and wrap input
				if (horizontal) {
					label.addClass('col-sm-3');
					var wrapper = angular.element('<div class="col-sm-9"></div>');
					wrapper.append(childElement);
				} else {

				}

				element.append(wrapper);

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
				if (passedName && !labelAdded) {
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