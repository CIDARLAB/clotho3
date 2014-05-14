'use strict';

describe('Controller: TestQuizCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestQuizCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestQuizCtrl = $controller('TestQuizCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
