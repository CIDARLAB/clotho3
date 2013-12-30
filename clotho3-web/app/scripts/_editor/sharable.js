angular.module('clotho.editor').controller('Editor_SharableCtrl', function($scope, $compile, Clotho) {

	$scope.$watch('sharable', function(newval) {
		console.log(newval);
		Clotho.get(newval.schema).then(function(result) {
			$scope.schema = result;
			$scope.schema_custom = result.custom;

			//todo - this is kinda a hack... should get element more directly
			var insert = angular.element(document).find('insert-fields').html(generateDynamicForm());
			$compile(insert.contents())($scope);
		});
	});


	/**
	 * @description Generates HTML for a form provided a properly formed schema. Implicit parameters are scope, with $scope.schema defined
	 * @returns {string}
	 */
	function generateDynamicForm () {
		var fulltext = "";

		//todo - need BSON types from server (GH #99)
		/*
		 string -> test vs. textarea
		 when use select?
		 what to do for objects?
		 id not editable
		 */
		var typeMap = {
			'string' : 'text',
			'boolean' : 'checkbox'
		};

		angular.forEach($scope.schema.fields, function(field) {


			var type = field.type || 'text';
			if (type == '?') field.type == 'text';
			var required = field.required ? "required='required'" : "";

			//todo - convert this to a directive (GH #105)
			var htmlText_pre = '<div class="control-group">' +
				'<label class="control-label" for="' + field.name + '">' + field.name + '</label>' +
				'<div class="controls">';
			var htmlText_post = '</div>' +
				'</div>';
			var inputText;

			switch (type) {
				case "textarea": {
					inputText = '<textarea class="input-large" id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-model="sharable.'+field.name+'" ng-disabled="!editMode"></textarea>';
					break;
				}
				case "select": {
					var optionsText = "";
					//todo - use ng-options
					angular.forEach(field.options, function(value, key) {
						optionsText = optionsText + '<option value="'+value+'">'+ value + '</option>';
					});

					inputText = '<select id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-disabled="!editMode" ng-model="sharable.'+field.name+'">' + optionsText + '</select>';
					break;
				}
				case "sharable": {
				}
				//todo - add filedrop support, and radio. checkbox works.
				default: {
					inputText = '<input type="' + type + '" class="input-large" id="' + field.name + '" name="' + field.name + '" ' + required + ' ng-disabled="!editMode" ng-model="sharable.'+field.name+'" >';
					break;
				}

			}

			fulltext += htmlText_pre + inputText + htmlText_post;
		});
		return fulltext;
	}


});