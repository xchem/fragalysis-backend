#!groovy​

// The 'xchem' Fragalysis Backend Jenkinsfile.

pipeline {

  agent { label 'buildah-slave' }

  environment {
    // Registry details
    USER = 'jenkins'
    REGISTRY = 'docker-registry.default:5000'
    STREAM_IMAGE = "${REGISTRY}/fragalysis-cicd/fragalysis-backend:latest"

    // Slack channel for all notifications
    SLACK_BUILD_CHANNEL = 'dls-builds'
    // Slack channel to be used for errors/failures
    SLACK_ALERT_CHANNEL = 'dls-alerts'
  }

  stages {

    stage('Inspect') {
      steps {
          slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                    message: "${JOB_NAME} build ${BUILD_NUMBER} - starting..."
          echo "Inspecting..."
          false
      }
    }

    stage('Build Image') {
      steps {
        echo "Building fragalysis-backend..."
        sh "buildah bud --format docker -f Dockerfile -t ${STREAM_IMAGE} ."
      }
    }

    stage('Push Image') {
      steps {
        script {
          TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
        }
        sh "podman login --tls-verify=false --username ${USER} --password ${TOKEN} ${REGISTRY}"
        sh "buildah push --tls-verify=false ${STREAM_IMAGE} docker://${STREAM_IMAGE}"
        sh "podman logout ${REGISTRY}"
      }
    }

  }

  // Post-job actions.
  // See https://jenkins.io/doc/book/pipeline/syntax/#post
  post {

    success {
      slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                color: 'good',
                message: "${JOB_NAME} build ${BUILD_NUMBER} - complete"
    }

    failure {
      slackSend channel: "#${SLACK_ALERT_CHANNEL}",
              color: 'danger',
              message: "${JOB_NAME} build ${BUILD_NUMBER} - failed (${BUILD_URL})"
    }

    fixed {
      slackSend channel: "#${SLACK_ALERT_CHANNEL}",
              color: 'good',
              message: "${JOB_NAME} build - fixed"
    }

  }

}
